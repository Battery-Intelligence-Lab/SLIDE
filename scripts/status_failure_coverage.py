"""Measure every literal ``Status`` failure arm in ``src/core``.

The gate deliberately works at failure-enumerator token granularity.  A line or
enclosing return statement can execute while the failure arm of a conditional
return remains unvisited, so ordinary line coverage is not sufficient.

Typical use in an instrumented build is::

    python scripts/status_failure_coverage.py prepare --build-dir build-coverage
    ctest --test-dir build-coverage --output-on-failure
    python scripts/status_failure_coverage.py report --build-dir build-coverage \
        --source-root . --exceptions tests/coverage/status_failure_exceptions.json
"""

from __future__ import annotations

import argparse
import bisect
import dataclasses
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Iterable, Sequence


FAILURE_CLASSES = {
    "defensive-only",
    "platform-specific",
    "unreachable-by-construction",
}
SOURCE_SUFFIXES = {".cc", ".cpp", ".cu", ".h", ".hh", ".hpp"}
SESSION_SCHEMA = 1
REPORT_SCHEMA = 1


class CoverageError(RuntimeError):
    """A deterministic coverage-gate failure."""


@dataclasses.dataclass(frozen=True)
class Token:
    value: str
    start: int
    end: int
    line: int
    column: int


@dataclasses.dataclass(frozen=True)
class Site:
    site_id: str
    path: str
    line: int
    column: int
    end_column: int
    status: str
    form: str
    statement: str
    source_sha256: str
    context_sha256: str


@dataclasses.dataclass(frozen=True)
class ManifestEntry:
    kind: str
    name: str
    binary: Path


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _line_starts(data: bytes) -> list[int]:
    starts = [0]
    starts.extend(index + 1 for index, value in enumerate(data) if value == 10)
    return starts


def _position(starts: Sequence[int], offset: int) -> tuple[int, int]:
    line_index = bisect.bisect_right(starts, offset) - 1
    return line_index + 1, offset - starts[line_index] + 1


def _skip_quoted(data: bytes, index: int, quote: int) -> int:
    index += 1
    while index < len(data):
        if data[index] == 92:  # backslash
            index += 2
        elif data[index] == quote:
            return index + 1
        else:
            index += 1
    raise CoverageError("unterminated C++ quoted literal")


def _skip_raw_string(data: bytes, index: int) -> int | None:
    prefix = next(
        (candidate for candidate in (b'u8R"', b'uR"', b'UR"', b'LR"', b'R"')
         if data.startswith(candidate, index)),
        None,
    )
    if prefix is None:
        return None
    # A raw-string delimiter is at most 16 bytes in C++.
    delimiter_start = index + len(prefix)
    open_paren = data.find(b"(", delimiter_start, min(len(data), delimiter_start + 17))
    if open_paren < 0:
        return None
    delimiter = data[delimiter_start:open_paren]
    terminator = b")" + delimiter + b'"'
    end = data.find(terminator, open_paren + 1)
    if end < 0:
        raise CoverageError("unterminated C++ raw string")
    return end + len(terminator)


def _tokens(data: bytes) -> list[Token]:
    """Tokenise enough C++ to locate return expressions without false text hits."""

    starts = _line_starts(data)
    tokens: list[Token] = []
    index = 0
    while index < len(data):
        value = data[index]
        if value in b" \t\r\n\f\v":
            index += 1
            continue
        if data[index : index + 2] == b"//":
            newline = data.find(b"\n", index + 2)
            index = len(data) if newline < 0 else newline + 1
            continue
        if data[index : index + 2] == b"/*":
            end = data.find(b"*/", index + 2)
            if end < 0:
                raise CoverageError("unterminated C++ block comment")
            index = end + 2
            continue
        raw_end = _skip_raw_string(data, index)
        if raw_end is not None:
            index = raw_end
            continue
        start = index
        if 48 <= value <= 57:
            # C++ digit separators use apostrophes (for example 131'072), so
            # consume a numeric preprocessing token before recognising chars.
            index += 1
            while index < len(data):
                nxt = data[index]
                if (
                    nxt in b"._'"
                    or nxt == 95
                    or 48 <= nxt <= 57
                    or 65 <= nxt <= 90
                    or 97 <= nxt <= 122
                ):
                    index += 1
                else:
                    break
        elif value in (34, 39):
            index = _skip_quoted(data, index, value)
            continue
        elif value == 95 or 65 <= value <= 90 or 97 <= value <= 122:
            index += 1
            while index < len(data):
                nxt = data[index]
                if nxt == 95 or 48 <= nxt <= 57 or 65 <= nxt <= 90 or 97 <= nxt <= 122:
                    index += 1
                else:
                    break
        elif data[index : index + 2] in (b"::", b"->", b"&&", b"||"):
            index += 2
        else:
            index += 1
        line, column = _position(starts, start)
        tokens.append(
            Token(data[start:index].decode("ascii", errors="replace"), start, index, line, column)
        )
    return tokens


def _matching_statement(tokens: Sequence[Token], return_index: int) -> tuple[int, list[Token]]:
    depths = {"(": 0, "[": 0, "{": 0}
    closing = {")": "(", "]": "[", "}": "{"}
    expression: list[Token] = []
    for index in range(return_index + 1, len(tokens)):
        token = tokens[index]
        if token.value == ";" and all(depth == 0 for depth in depths.values()):
            return index, expression
        if token.value in depths:
            depths[token.value] += 1
        elif token.value in closing:
            opener = closing[token.value]
            if depths[opener] == 0:
                raise CoverageError(
                    f"unbalanced return expression at {tokens[return_index].line}:"
                    f"{tokens[return_index].column}"
                )
            depths[opener] -= 1
        expression.append(token)
    raise CoverageError(
        f"unterminated return statement at {tokens[return_index].line}:"
        f"{tokens[return_index].column}"
    )


def _status_sequences(tokens: Sequence[Token]) -> list[tuple[int, Token]]:
    result: list[tuple[int, Token]] = []
    for index in range(len(tokens) - 2):
        if (
            tokens[index].value == "Status"
            and tokens[index + 1].value == "::"
            and re.fullmatch(r"[A-Za-z_]\w*", tokens[index + 2].value)
        ):
            result.append((index, tokens[index + 2]))
    return result


def _strip_outer_parentheses(tokens: Sequence[Token]) -> list[Token]:
    result = list(tokens)
    while len(result) >= 2 and result[0].value == "(" and result[-1].value == ")":
        depth = 0
        closes_at_end = False
        for index, token in enumerate(result):
            if token.value == "(":
                depth += 1
            elif token.value == ")":
                depth -= 1
                if depth == 0:
                    closes_at_end = index == len(result) - 1
                    break
        if not closes_at_end:
            break
        result = result[1:-1]
    return result


def _status_literal(tokens: Sequence[Token]) -> Token | None:
    stripped = _strip_outer_parentheses(tokens)
    values = [token.value for token in stripped]
    if len(values) == 3 and values[0] == "Status" and values[1] == "::":
        return stripped[2]
    if (
        len(values) == 5
        and values[0] == "slide"
        and values[1] == "::"
        and values[2] == "Status"
        and values[3] == "::"
    ):
        return stripped[4]
    return None


def _conditional_arms(tokens: Sequence[Token]) -> tuple[list[Token], list[Token]] | None:
    tokens = _strip_outer_parentheses(tokens)
    round_depth = square_depth = brace_depth = 0
    question_index: int | None = None
    colon_index: int | None = None
    for index, token in enumerate(tokens):
        if token.value == "(":
            round_depth += 1
        elif token.value == ")":
            round_depth -= 1
        elif token.value == "[":
            square_depth += 1
        elif token.value == "]":
            square_depth -= 1
        elif token.value == "{":
            brace_depth += 1
        elif token.value == "}":
            brace_depth -= 1
        elif (
            token.value == "?"
            and round_depth == 0
            and square_depth == 0
            and brace_depth == 0
        ):
            if question_index is not None:
                raise CoverageError("nested or chained top-level conditional return is unsupported")
            question_index = index
        elif (
            token.value == ":"
            and question_index is not None
            and round_depth == 0
            and square_depth == 0
            and brace_depth == 0
        ):
            if colon_index is not None:
                raise CoverageError("conditional return has multiple top-level colons")
            colon_index = index
    if question_index is None:
        return None
    if colon_index is None or colon_index < question_index:
        raise CoverageError("conditional return has no matching top-level colon")
    return tokens[question_index + 1 : colon_index], tokens[colon_index + 1 :]


def _normalise_statement(data: bytes, start: int, end: int) -> str:
    return " ".join(data[start:end].decode("utf-8", errors="replace").split())


def scan_sites(source_root: Path) -> list[Site]:
    core = source_root / "src" / "core"
    sites: list[Site] = []
    seen_ids: set[str] = set()
    for path in sorted(item for item in core.rglob("*") if item.suffix in SOURCE_SUFFIXES):
        data = path.read_bytes()
        try:
            tokens = _tokens(data)
        except CoverageError as error:
            raise CoverageError(f"{path}: {error}") from error
        starts = _line_starts(data)
        lines = data.splitlines(keepends=True)
        for return_index, token in enumerate(tokens):
            if token.value != "return":
                continue
            semicolon_index, expression = _matching_statement(tokens, return_index)
            del semicolon_index  # the token index is not otherwise needed
            sequences = _status_sequences(expression)
            failures = [(index, enum) for index, enum in sequences if enum.value != "Success"]
            if not failures:
                continue
            if any(item.value == "return" for item in expression):
                raise CoverageError(
                    f"failure literal inside a nested return at {path}:{token.line}:"
                    f"{token.column} is ambiguous"
                )
            arms = _conditional_arms(expression)
            conditional = arms is not None
            if not conditional:
                literal = _status_literal(expression)
                if literal is None or literal.value != failures[0][1].value or len(failures) != 1:
                    raise CoverageError(
                        f"non-literal direct Status return at {path}:{token.line}:{token.column}"
                    )
            else:
                assert arms is not None
                arm_literals = [candidate for arm in arms if (candidate := _status_literal(arm))]
                failure_tokens = {enum.start for _, enum in failures}
                literal_failures = {
                    enum.start for enum in arm_literals if enum.value != "Success"
                }
                if literal_failures != failure_tokens:
                    raise CoverageError(
                        f"non-literal conditional Status arm at {path}:{token.line}:"
                        f"{token.column}"
                    )
            statement_end = expression[-1].end if expression else token.end
            statement = _normalise_statement(data, token.start, statement_end)
            statement_hash = _sha256(data[token.start:statement_end])
            first_context_line = max(0, token.line - 3)
            last_line = max(enum.line for _, enum in failures)
            context = b"".join(lines[first_context_line : min(len(lines), last_line + 2)])
            context_hash = _sha256(context)
            relative = path.relative_to(source_root).as_posix()
            for _, enum in failures:
                site_id = f"{relative}:{enum.line}:{enum.column}:{enum.value}"
                if site_id in seen_ids:
                    raise CoverageError(f"duplicate failure site ID: {site_id}")
                seen_ids.add(site_id)
                sites.append(
                    Site(
                        site_id=site_id,
                        path=relative,
                        line=enum.line,
                        column=enum.column,
                        end_column=enum.column + len(enum.value.encode("ascii")),
                        status=enum.value,
                        form="conditional" if conditional else "direct",
                        statement=statement,
                        source_sha256=statement_hash,
                        context_sha256=context_hash,
                    )
                )
    return sorted(sites, key=lambda site: (site.path, site.line, site.column, site.status))


def read_manifest(build_dir: Path) -> list[ManifestEntry]:
    path = build_dir / "coverage" / "targets.tsv"
    if not path.is_file():
        raise CoverageError(f"coverage target manifest is missing: {path}")
    entries: list[ManifestEntry] = []
    seen_names: set[str] = set()
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 3 or fields[0] not in {"anchor", "test"}:
            raise CoverageError(f"invalid target manifest line {number}: {line!r}")
        kind, name, binary_text = fields
        binary = Path(binary_text).resolve()
        if name in seen_names:
            raise CoverageError(f"duplicate coverage target name: {name}")
        if not binary.is_file():
            raise CoverageError(f"coverage binary does not exist: {binary}")
        seen_names.add(name)
        entries.append(ManifestEntry(kind, name, binary))
    anchors = [entry for entry in entries if entry.kind == "anchor"]
    tests = [entry for entry in entries if entry.kind == "test"]
    if len(anchors) != 1 or not tests:
        raise CoverageError(
            "coverage manifest must contain exactly one anchor and at least one test"
        )
    return entries


def _safe_child(path: Path, parent: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(parent.resolve())
    except ValueError as error:
        raise CoverageError(f"refusing to modify path outside build tree: {resolved}") from error
    return resolved


def _source_identities(source_root: Path) -> dict[str, str]:
    core = source_root / "src" / "core"
    return {
        path.relative_to(source_root).as_posix(): _sha256(path.read_bytes())
        for path in sorted(item for item in core.rglob("*") if item.suffix in SOURCE_SUFFIXES)
    }


def _compiled_source_identities(build_dir: Path, source_root: Path) -> dict[str, str]:
    entries = json.loads((build_dir / "compile_commands.json").read_text(encoding="utf-8"))
    roots = (
        source_root / "src",
        source_root / "tests" / "unit",
        source_root / "tests" / "parity",
        source_root / "tests" / "coverage",
    )
    result: dict[str, str] = {}
    for entry in entries:
        path = Path(str(entry["file"])).resolve()
        if not any(path.is_relative_to(root.resolve()) for root in roots):
            continue
        if not path.is_file():
            raise CoverageError(f"compiled source does not exist: {path}")
        result[path.relative_to(source_root).as_posix()] = _sha256(path.read_bytes())
    return dict(sorted(result.items()))


def _cache_value(build_dir: Path, name: str) -> str:
    prefix = name + ":"
    for line in (build_dir / "CMakeCache.txt").read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix) and "=" in line:
            return line.split("=", 1)[1]
    raise CoverageError(f"CMake cache value is missing: {name}")


def _coverage_tool(build_dir: Path, override: str | None, cache_name: str) -> str:
    """Use an explicit tool or the exact version selected at CMake configure time."""
    return override if override is not None else _cache_value(build_dir, cache_name)


def _require_up_to_date(build_dir: Path) -> None:
    ninja = _cache_value(build_dir, "CMAKE_MAKE_PROGRAM")
    result = subprocess.run(
        [ninja, "-C", str(build_dir), "-n", "slide_coverage_binaries"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    output = (result.stdout + result.stderr).decode("utf-8", errors="replace").strip()
    if result.returncode != 0:
        raise CoverageError(f"coverage freshness probe failed: {output}")
    if "ninja: no work to do." not in output:
        raise CoverageError(
            "coverage binaries are out of date; build slide_coverage_binaries first\n"
            + output
        )


def prepare(build_dir: Path, source_root: Path) -> None:
    build_dir = build_dir.resolve()
    source_root = source_root.resolve()
    entries = read_manifest(build_dir)
    _require_up_to_date(build_dir)
    coverage_dir = _safe_child(build_dir / "coverage", build_dir)
    raw_dir = _safe_child(coverage_dir / "raw", build_dir)
    work_dir = _safe_child(coverage_dir / "work", build_dir)
    shutil.rmtree(raw_dir, ignore_errors=True)
    shutil.rmtree(work_dir, ignore_errors=True)
    raw_dir.mkdir(parents=True)
    work_dir.mkdir(parents=True)
    started_ns = time.time_ns()
    session = {
        "schema_version": SESSION_SCHEMA,
        "started_ns": started_ns,
        "manifest_sha256": _sha256((coverage_dir / "targets.tsv").read_bytes()),
        "compile_commands_sha256": _sha256((build_dir / "compile_commands.json").read_bytes()),
        "compiled_sources": _compiled_source_identities(build_dir, source_root),
        "source_root": str(source_root),
        "sources": _source_identities(source_root),
        "binaries": {
            entry.name: {
                "size": entry.binary.stat().st_size,
                "mtime_ns": entry.binary.stat().st_mtime_ns,
            }
            for entry in entries
        },
    }
    (coverage_dir / "session.json").write_text(
        json.dumps(session, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"prepared fresh coverage session for {len(entries) - 1} tests")


def _load_session(
    build_dir: Path, source_root: Path, entries: Sequence[ManifestEntry]
) -> dict[str, object]:
    coverage_dir = build_dir / "coverage"
    session_path = coverage_dir / "session.json"
    if not session_path.is_file():
        raise CoverageError("coverage session is missing; run the prepare command first")
    session = json.loads(session_path.read_text(encoding="utf-8"))
    if session.get("schema_version") != SESSION_SCHEMA:
        raise CoverageError("unsupported coverage session schema")
    if session.get("manifest_sha256") != _sha256((coverage_dir / "targets.tsv").read_bytes()):
        raise CoverageError("coverage manifest changed after the session was prepared")
    if session.get("compile_commands_sha256") != _sha256(
        (build_dir / "compile_commands.json").read_bytes()
    ):
        raise CoverageError("compile_commands.json changed after prepare")
    if session.get("source_root") != str(source_root.resolve()):
        raise CoverageError("coverage source root differs from the prepared session")
    if session.get("sources") != _source_identities(source_root.resolve()):
        raise CoverageError("src/core changed after the coverage session was prepared")
    if session.get("compiled_sources") != _compiled_source_identities(
        build_dir, source_root.resolve()
    ):
        raise CoverageError("compiled source changed after the coverage session was prepared")
    recorded = session.get("binaries")
    if not isinstance(recorded, dict):
        raise CoverageError("coverage session has no binary identity records")
    for entry in entries:
        identity = recorded.get(entry.name)
        if not isinstance(identity, dict):
            raise CoverageError(f"coverage session omitted binary {entry.name}")
        stat = entry.binary.stat()
        if identity.get("size") != stat.st_size or identity.get("mtime_ns") != stat.st_mtime_ns:
            raise CoverageError(f"coverage binary changed after prepare: {entry.name}")
    return session


def _profile_groups(
    build_dir: Path, entries: Sequence[ManifestEntry], started_ns: int
) -> dict[str, list[Path]]:
    tests = [entry for entry in entries if entry.kind == "test"]
    groups = {entry.name: [] for entry in tests}
    raw_files = sorted((build_dir / "coverage" / "raw").glob("*.profraw"))
    if not raw_files:
        raise CoverageError("no raw profiles were produced")
    ordered_names = sorted(groups, key=len, reverse=True)
    for profile in raw_files:
        matches = [name for name in ordered_names if profile.name.startswith(name + "-")]
        if len(matches) != 1:
            raise CoverageError(f"unexpected or ambiguous raw profile: {profile.name}")
        # DrvFS and some hosted filesystems expose coarse timestamps.  Cleaning
        # the directory is the primary freshness guarantee; allow one second.
        if profile.stat().st_mtime_ns + 1_000_000_000 < started_ns:
            raise CoverageError(f"stale raw profile: {profile.name}")
        groups[matches[0]].append(profile)
    missing = [name for name, profiles in groups.items() if not profiles]
    if missing:
        raise CoverageError("tests produced no profile: " + ", ".join(sorted(missing)))
    return groups


def _run(command: Sequence[str], *, cwd: Path, stdout: Path | None = None) -> str:
    output_handle = stdout.open("wb") if stdout is not None else subprocess.PIPE
    try:
        result = subprocess.run(
            list(command),
            cwd=cwd,
            stdout=output_handle,
            stderr=subprocess.PIPE,
            check=False,
        )
    finally:
        if stdout is not None:
            output_handle.close()
    stderr = result.stderr.decode("utf-8", errors="replace").strip()
    if result.returncode != 0:
        raise CoverageError(
            f"command failed ({result.returncode}): {shlex.join(command)}"
            + (f"\n{stderr}" if stderr else "")
        )
    if stderr:
        raise CoverageError(f"coverage tool emitted diagnostics: {stderr}")
    if stdout is None:
        assert isinstance(result.stdout, bytes)
        return result.stdout.decode("utf-8", errors="replace")
    return ""


def _tool_major(tool: str, cwd: Path) -> tuple[int, str]:
    text = _run([tool, "--version"], cwd=cwd)
    match = re.search(r"(?:version\s+)(\d+)", text)
    if not match:
        raise CoverageError(f"cannot parse tool version from {tool}: {text.splitlines()[:1]}")
    return int(match.group(1)), text.splitlines()[0].strip()


def _merge_profiles(tool: str, profiles: Iterable[Path], output: Path, cwd: Path) -> None:
    profile_list = [str(path) for path in profiles]
    if not profile_list:
        raise CoverageError("cannot merge an empty profile set")
    _run([tool, "merge", "-sparse", *profile_list, "-o", str(output)], cwd=cwd)


def _run_anchor(anchor: Path, work_dir: Path, cwd: Path) -> list[Path]:
    for stale in work_dir.glob("anchor-*.profraw"):
        stale.unlink()
    pattern = work_dir / "anchor-%p-%m.profraw"
    environment = os.environ.copy()
    environment["LLVM_PROFILE_FILE"] = str(pattern)
    result = subprocess.run(
        [str(anchor)],
        cwd=cwd,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    diagnostics = (result.stdout + result.stderr).decode("utf-8", errors="replace").strip()
    if result.returncode != 0 or diagnostics:
        raise CoverageError(
            f"coverage anchor failed ({result.returncode})"
            + (f": {diagnostics}" if diagnostics else "")
        )
    profiles = sorted(work_dir.glob("anchor-*.profraw"))
    if not profiles:
        raise CoverageError("coverage anchor produced no zero-count profile")
    return profiles


def _export(
    llvm_cov: str,
    objects: Sequence[Path],
    profile: Path,
    output: Path,
    cwd: Path,
) -> None:
    if not objects:
        raise CoverageError("coverage export has no object")
    command = [
        llvm_cov,
        "export",
        str(objects[0]),
        f"-instr-profile={profile}",
        "--skip-functions",
        "--skip-expansions",
    ]
    for path in objects[1:]:
        command.append(f"-object={path}")
    _run(command, cwd=cwd, stdout=output)


def _coverage_files(export_path: Path) -> dict[str, list[list[list[object]]]]:
    payload = json.loads(export_path.read_text(encoding="utf-8"))
    if payload.get("type") != "llvm.coverage.json.export":
        raise CoverageError("unexpected llvm-cov JSON type")
    version = str(payload.get("version", ""))
    if not version.startswith("2."):
        raise CoverageError(f"unsupported llvm-cov JSON schema: {version}")
    result: dict[str, list[list[list[object]]]] = {}
    for data in payload.get("data", []):
        for item in data.get("files", []):
            filename = str(Path(item["filename"]).resolve())
            result.setdefault(filename, []).append(item.get("segments", []))
    return result


def _segment_interval_count(
    segments: Sequence[Sequence[object]], line: int, column: int, end_column: int
) -> int | None:
    grouped: dict[tuple[int, int], Sequence[object]] = {}
    for segment in segments:
        if len(segment) < 6:
            raise CoverageError("malformed llvm-cov segment")
        grouped[(int(segment[0]), int(segment[1]))] = segment
    coordinates = sorted(grouped)
    start = (line, column)
    index = bisect.bisect_right(coordinates, start) - 1
    if index < 0:
        return None
    coordinate = coordinates[index]
    segment = grouped[coordinate]
    if not bool(segment[3]) or bool(segment[5]):
        return None
    next_coordinate = coordinates[index + 1] if index + 1 < len(coordinates) else None
    token_end = (line, end_column)
    if next_coordinate is not None and next_coordinate < token_end:
        raise CoverageError(
            f"failure token crosses coverage segment boundary at {line}:{column}"
        )
    return int(segment[2])


def _counts_for_sites(
    files: dict[str, list[list[list[object]]]], sites: Sequence[Site], source_root: Path
) -> dict[str, int | None]:
    counts: dict[str, int | None] = {}
    for site in sites:
        filename = str((source_root / site.path).resolve())
        mappings = files.get(filename, [])
        observed: list[int] = []
        for segments in mappings:
            count = _segment_interval_count(segments, site.line, site.column, site.end_column)
            if count is not None:
                observed.append(count)
        counts[site.site_id] = sum(observed) if observed else None
    return counts


def _compile_arguments(entry: dict[str, object]) -> list[str]:
    arguments = entry.get("arguments")
    if isinstance(arguments, list) and all(isinstance(item, str) for item in arguments):
        return list(arguments)
    command = entry.get("command")
    if isinstance(command, str):
        return shlex.split(command)
    raise CoverageError("compile_commands entry has neither arguments nor command")


def _preprocess_command(arguments: Sequence[str], original: Path, target: Path) -> list[str]:
    result = [arguments[0]]
    skip_next = False
    options_with_value = {"-o", "-MF", "-MT", "-MQ", "--serialize-diagnostics"}
    for argument in arguments[1:]:
        if skip_next:
            skip_next = False
            continue
        if argument in options_with_value:
            skip_next = True
            continue
        if argument in {"-c", "-MD", "-MMD", "-MP"}:
            continue
        try:
            if Path(argument).resolve() == original.resolve():
                continue
        except (OSError, ValueError):
            pass
        result.append(argument)
    result.extend(["-E", "-fdirectives-only", "-x", "c++", str(target)])
    return result


def _active_lines_for_file(
    path: Path, compile_entry: dict[str, object], original: Path
) -> set[int]:
    arguments = _compile_arguments(compile_entry)
    command = _preprocess_command(arguments, original, path)
    directory = Path(str(compile_entry["directory"]))
    process = subprocess.Popen(
        command,
        cwd=directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert process.stdout is not None
    current_file: Path | None = None
    current_line = 0
    active: set[int] = set()
    marker = re.compile(r'^\s*#\s+(\d+)\s+"([^"]+)"')
    target = path.resolve()
    for output_line in process.stdout:
        match = marker.match(output_line)
        if match:
            current_line = int(match.group(1))
            current_file = Path(match.group(2)).resolve()
            continue
        if (
            current_file == target
            and output_line.strip()
            and not output_line.lstrip().startswith("#")
        ):
            active.add(current_line)
        current_line += 1
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise CoverageError(
            f"preprocessing failed for {path} ({return_code}):\n{stderr.strip()}"
        )
    return active


def classify_active_sites(
    sites: Sequence[Site], source_root: Path, build_dir: Path
) -> tuple[dict[str, bool], str]:
    compile_commands_path = build_dir / "compile_commands.json"
    if not compile_commands_path.is_file():
        raise CoverageError("compile_commands.json is required for activity classification")
    entries = json.loads(compile_commands_path.read_text(encoding="utf-8"))
    if not isinstance(entries, list) or not entries:
        raise CoverageError("compile_commands.json is empty")
    by_file: dict[Path, dict[str, object]] = {}
    for entry in entries:
        if not isinstance(entry, dict) or "file" not in entry:
            raise CoverageError("malformed compile_commands.json entry")
        by_file[Path(str(entry["file"])).resolve()] = entry
    base_original, base_entry = next(
        (
            (path, entry)
            for path, entry in by_file.items()
            if "src/core/" in path.as_posix() and path.suffix in {".cc", ".cpp"}
        ),
        (None, None),
    )
    if base_entry is None or base_original is None:
        raise CoverageError("compile database contains no core production translation unit")
    compiler = _compile_arguments(base_entry)[0]
    active: dict[str, bool] = {}
    paths = sorted({site.path for site in sites})
    for relative in paths:
        path = (source_root / relative).resolve()
        entry = by_file.get(path, base_entry)
        original = path if path in by_file else base_original
        active_lines = _active_lines_for_file(path, entry, original)
        for site in (candidate for candidate in sites if candidate.path == relative):
            active[site.site_id] = site.line in active_lines
    return active, compiler


def _load_exceptions(path: Path, sites: Sequence[Site]) -> dict[str, dict[str, str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or not isinstance(payload.get("exceptions"), list):
        raise CoverageError("invalid exception-manifest schema")
    entries = payload["exceptions"]
    if len(entries) > 10:
        raise CoverageError(f"exception budget exceeded: {len(entries)} > 10")
    known = {site.site_id: site for site in sites}
    result: dict[str, dict[str, str]] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise CoverageError("exception entry must be an object")
        required = {"site_id", "source_sha256", "context_sha256", "class", "reason"}
        if set(entry) != required:
            raise CoverageError(
                "exception entry must contain exactly " + ", ".join(sorted(required))
            )
        site_id = str(entry["site_id"])
        if any(character in site_id for character in "*?[]{}") or ".." in site_id:
            raise CoverageError(f"wildcard or range exception is forbidden: {site_id}")
        if site_id in result:
            raise CoverageError(f"duplicate exception: {site_id}")
        site = known.get(site_id)
        if site is None:
            raise CoverageError(f"stale exception site: {site_id}")
        if entry["source_sha256"] != site.source_sha256:
            raise CoverageError(f"stale source hash for exception: {site_id}")
        if entry["context_sha256"] != site.context_sha256:
            raise CoverageError(f"stale context hash for exception: {site_id}")
        if entry["class"] not in FAILURE_CLASSES:
            raise CoverageError(f"invalid exception class for {site_id}: {entry['class']}")
        reason = str(entry["reason"]).strip()
        if len(reason) < 40:
            raise CoverageError(f"exception reason is not substantive: {site_id}")
        result[site_id] = {
            "class": str(entry["class"]),
            "reason": reason,
        }
    return result


def _write_markdown(path: Path, report: dict[str, object]) -> None:
    totals = report["totals"]
    assert isinstance(totals, dict)
    sites = report["sites"]
    assert isinstance(sites, list)
    per_file: dict[str, dict[str, int]] = {}
    for site in sites:
        assert isinstance(site, dict)
        summary = per_file.setdefault(
            str(site["path"]), {"active": 0, "covered": 0, "excepted": 0, "inactive": 0}
        )
        disposition = str(site["disposition"])
        if disposition == "inactive":
            summary["inactive"] += 1
        else:
            summary["active"] += 1
            if disposition == "covered":
                summary["covered"] += 1
            elif disposition == "excepted":
                summary["excepted"] += 1
    lines = [
        "# P9-G3 Status-failure coverage",
        "",
        f"Generated: {report['generated_utc']}",
        "",
        f"Result: **{report['result']}**",
        "",
        f"- Lexical failure arms: {totals['lexical']}",
        f"- Active optional-off arms: {totals['active']}",
        f"- Measured covered arms: {totals['covered']}",
        f"- Structural exceptions: {totals['excepted']} / 10",
        f"- Inactive optional branches: {totals['inactive']}",
        f"- Uncovered/unmapped active arms: {totals['failed']}",
        "",
        "| File | Active | Covered | Excepted | Inactive |",
        "|---|---:|---:|---:|---:|",
    ]
    for filename, summary in sorted(per_file.items()):
        lines.append(
            f"| `{filename}` | {summary['active']} | {summary['covered']} | "
            f"{summary['excepted']} | {summary['inactive']} |"
        )
    exceptions = [site for site in sites if site["disposition"] == "excepted"]
    failures = [site for site in sites if site["disposition"] in {"uncovered", "unmapped"}]
    lines.extend(["", "## Structural exceptions", ""])
    if exceptions:
        for site in exceptions:
            lines.append(
                f"- `{site['site_id']}` — **{site['exception_class']}**: {site['reason']}"
            )
    else:
        lines.append("None.")
    lines.extend(["", "## Gate failures", ""])
    if failures:
        for site in failures:
            lines.append(f"- `{site['site_id']}` — {site['disposition']}")
    else:
        lines.append("None.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def report(
    source_root: Path,
    build_dir: Path,
    exceptions_path: Path,
    json_output: Path,
    markdown_output: Path,
    llvm_cov: str,
    llvm_profdata: str,
) -> None:
    source_root = source_root.resolve()
    build_dir = build_dir.resolve()
    entries = read_manifest(build_dir)
    session = _load_session(build_dir, source_root, entries)
    started_ns = int(session["started_ns"])
    profile_groups = _profile_groups(build_dir, entries, started_ns)
    sites = scan_sites(source_root)
    active, compiler = classify_active_sites(sites, source_root, build_dir)
    compiler_major, compiler_version = _tool_major(compiler, build_dir)
    cov_major, cov_version = _tool_major(llvm_cov, build_dir)
    profdata_major, profdata_version = _tool_major(llvm_profdata, build_dir)
    if len({compiler_major, cov_major, profdata_major}) != 1:
        raise CoverageError(
            "coverage tool major versions differ: "
            f"compiler={compiler_major}, llvm-cov={cov_major}, llvm-profdata={profdata_major}"
        )
    exceptions = _load_exceptions(exceptions_path, sites)
    work_dir = build_dir / "coverage" / "work"
    anchor = next(entry for entry in entries if entry.kind == "anchor")
    tests = [entry for entry in entries if entry.kind == "test"]
    # Different standalone test binaries can contain same-name inline/template
    # functions with different mapping hashes. Export each test only against
    # its own profile. A separately executed whole-archive anchor supplies the
    # zero-count mapping union for production .cpp members that no test links.
    anchor_profile = work_dir / "anchor.profdata"
    anchor_export = work_dir / "anchor.json"
    _merge_profiles(
        llvm_profdata,
        _run_anchor(anchor.binary, work_dir, build_dir),
        anchor_profile,
        build_dir,
    )
    _export(llvm_cov, [anchor.binary], anchor_profile, anchor_export, build_dir)
    anchor_counts = _counts_for_sites(_coverage_files(anchor_export), sites, source_root)
    per_test: dict[str, dict[str, int | None]] = {}
    for entry in tests:
        profile = work_dir / f"{entry.name}.profdata"
        export_path = work_dir / f"{entry.name}.json"
        _merge_profiles(llvm_profdata, profile_groups[entry.name], profile, build_dir)
        _export(llvm_cov, [entry.binary], profile, export_path, build_dir)
        per_test[entry.name] = _counts_for_sites(
            _coverage_files(export_path), sites, source_root
        )

    site_records: list[dict[str, object]] = []
    errors: list[str] = []
    for site in sites:
        is_active = active[site.site_id]
        hit_by = {
            name: int(count)
            for name, counts in per_test.items()
            if (count := counts[site.site_id]) is not None and count > 0
        }
        observed_count = sum(hit_by.values())
        exception = exceptions.get(site.site_id)
        mapped = anchor_counts[site.site_id] is not None or any(
            counts[site.site_id] is not None for counts in per_test.values()
        )
        if not is_active:
            disposition = "inactive"
            if exception is not None:
                errors.append(f"inactive site has a stale exception: {site.site_id}")
        elif not mapped:
            disposition = "unmapped"
            errors.append(f"active site has no coverage mapping: {site.site_id}")
        elif observed_count > 0:
            disposition = "covered"
            if exception is not None:
                errors.append(f"covered site has a stale exception: {site.site_id}")
        elif exception is not None:
            disposition = "excepted"
        else:
            disposition = "uncovered"
            errors.append(f"active failure arm is uncovered: {site.site_id}")
        site_records.append(
            {
                "site_id": site.site_id,
                "path": site.path,
                "line": site.line,
                "column": site.column,
                "status": site.status,
                "form": site.form,
                "active": is_active,
                "mapped": mapped,
                "observed_count": observed_count,
                "hit_by": hit_by,
                "disposition": disposition,
                "exception_class": None if exception is None else exception["class"],
                "reason": None if exception is None else exception["reason"],
                "source_sha256": site.source_sha256,
                "context_sha256": site.context_sha256,
                "statement": site.statement,
            }
        )
    dispositions = [str(record["disposition"]) for record in site_records]
    totals = {
        "lexical": len(site_records),
        "active": sum(bool(record["active"]) for record in site_records),
        "covered": dispositions.count("covered"),
        "excepted": dispositions.count("excepted"),
        "inactive": dispositions.count("inactive"),
        "failed": dispositions.count("uncovered") + dispositions.count("unmapped"),
    }
    payload: dict[str, object] = {
        "schema_version": REPORT_SCHEMA,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "result": "PASS" if not errors else "FAIL",
        "configuration": "Linux Clang optional-off",
        "tools": {
            "compiler": compiler_version,
            "llvm_cov": cov_version,
            "llvm_profdata": profdata_version,
        },
        "totals": totals,
        "errors": errors,
        "sites": site_records,
    }
    json_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_markdown(markdown_output, payload)
    print(
        f"Status failure coverage: {payload['result']} — {totals['covered']} covered, "
        f"{totals['excepted']} excepted, {totals['failed']} failed, "
        f"{totals['inactive']} inactive"
    )
    if errors:
        raise CoverageError("\n".join(errors))


def self_test() -> None:
    sample = b'''// return slide::Status::Invalid_states;\n
const char *s = "return slide::Status::Numerical_failure;";
const char *r = R"tag(return slide::Status::Invalid_parameters;)tag";
const char8_t *r8 = u8R"tag(" return slide::Status::Invalid_states;)tag";
slide::Status f(bool ok) {
  if (!ok) return slide::Status::Invalid_parameters;
  return ok ? slide::Status::Success
            : slide::Status::Numerical_failure;
}
'''
    tokens = _tokens(sample)
    returns = [index for index, token in enumerate(tokens) if token.value == "return"]
    if len(returns) != 2:
        raise CoverageError(f"lexer self-test found {len(returns)} real returns, expected 2")
    expressions = [_matching_statement(tokens, index)[1] for index in returns]
    failures = [
        enum.value
        for expression in expressions
        for _, enum in _status_sequences(expression)
        if enum.value != "Success"
    ]
    if failures != ["Invalid_parameters", "Numerical_failure"]:
        raise CoverageError(f"lexer self-test found wrong failures: {failures}")
    segments = [
        [1, 1, 3, True, True, False],
        [1, 8, 0, True, True, False],
        [1, 16, 0, False, False, False],
    ]
    if _segment_interval_count(segments, 1, 2, 7) != 3:
        raise CoverageError("segment self-test missed covered interval")
    if _segment_interval_count(segments, 1, 9, 15) != 0:
        raise CoverageError("segment self-test confused zero count with unmapped")
    if _segment_interval_count(segments, 1, 17, 20) is not None:
        raise CoverageError("segment self-test mapped an uncovered gap")
    try:
        _segment_interval_count(segments, 1, 7, 9)
    except CoverageError:
        pass
    else:
        raise CoverageError("segment self-test accepted a token crossing a boundary")

    utf8 = "void f() { /* μ */ return slide::Status::Invalid_states; }\r\n".encode()
    utf8_tokens = _tokens(utf8)
    enum = next(token for token in utf8_tokens if token.value == "Invalid_states")
    expected_column = utf8.index(b"Invalid_states") + 1
    if enum.line != 1 or enum.column != expected_column:
        raise CoverageError("lexer self-test does not use UTF-8 byte columns")

    fake_site = Site(
        site_id="src/core/Fake.cpp:1:1:Invalid_states",
        path="src/core/Fake.cpp",
        line=1,
        column=1,
        end_column=15,
        status="Invalid_states",
        form="direct",
        statement="return slide::Status::Invalid_states",
        source_sha256="a" * 64,
        context_sha256="b" * 64,
    )
    valid_exception = {
        "site_id": fake_site.site_id,
        "source_sha256": fake_site.source_sha256,
        "context_sha256": fake_site.context_sha256,
        "class": "defensive-only",
        "reason": "A deterministic self-test reason that is deliberately substantive.",
    }
    with tempfile.TemporaryDirectory() as directory_text:
        directory = Path(directory_text)
        (directory / "CMakeCache.txt").write_text(
            "SLIDE_LLVM_COV_EXECUTABLE:FILEPATH=/opt/llvm-cov-18\n"
            "SLIDE_LLVM_PROFDATA_EXECUTABLE:FILEPATH=/opt/llvm-profdata-18\n",
            encoding="utf-8",
        )
        if _coverage_tool(
            directory, None, "SLIDE_LLVM_COV_EXECUTABLE"
        ) != "/opt/llvm-cov-18":
            raise CoverageError("tool self-test ignored the configured llvm-cov")
        if _coverage_tool(
            directory, "/explicit/llvm-profdata", "SLIDE_LLVM_PROFDATA_EXECUTABLE"
        ) != "/explicit/llvm-profdata":
            raise CoverageError("tool self-test ignored an explicit override")
        manifest = directory / "exceptions.json"
        manifest.write_text(
            json.dumps({"schema_version": 1, "exceptions": [valid_exception]}),
            encoding="utf-8",
        )
        if fake_site.site_id not in _load_exceptions(manifest, [fake_site]):
            raise CoverageError("exception self-test rejected an exact entry")
        wildcard = dict(valid_exception, site_id="src/core/Fake.cpp:*")
        manifest.write_text(
            json.dumps({"schema_version": 1, "exceptions": [wildcard]}),
            encoding="utf-8",
        )
        try:
            _load_exceptions(manifest, [fake_site])
        except CoverageError:
            pass
        else:
            raise CoverageError("exception self-test accepted a wildcard")
        manifest.write_text(
            json.dumps({"schema_version": 1, "exceptions": [valid_exception] * 11}),
            encoding="utf-8",
        )
        try:
            _load_exceptions(manifest, [fake_site])
        except CoverageError:
            pass
        else:
            raise CoverageError("exception self-test exceeded the ten-site budget")
    print("status-failure coverage self-test: PASS")


def census(source_root: Path, build_dir: Path | None) -> None:
    source_root = source_root.resolve()
    sites = scan_sites(source_root)
    direct = sum(site.form == "direct" for site in sites)
    conditional = len(sites) - direct
    statuses: dict[str, int] = {}
    for site in sites:
        statuses[site.status] = statuses.get(site.status, 0) + 1
    payload: dict[str, object] = {
        "conditional": conditional,
        "direct": direct,
        "statuses": dict(sorted(statuses.items())),
        "total": len(sites),
    }
    if build_dir is not None:
        activity, compiler = classify_active_sites(sites, source_root, build_dir.resolve())
        payload.update(
            {
                "active": sum(activity.values()),
                "compiler": compiler,
                "inactive": len(activity) - sum(activity.values()),
            }
        )
    print(json.dumps(payload, sort_keys=True))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare", help="start a fresh profile session")
    prepare_parser.add_argument("--build-dir", type=Path, required=True)
    prepare_parser.add_argument("--source-root", type=Path, default=Path.cwd())
    report_parser = subparsers.add_parser("report", help="measure and enforce the gate")
    report_parser.add_argument("--source-root", type=Path, default=Path.cwd())
    report_parser.add_argument("--build-dir", type=Path, required=True)
    report_parser.add_argument("--exceptions", type=Path, required=True)
    report_parser.add_argument("--json-output", type=Path, required=True)
    report_parser.add_argument("--markdown-output", type=Path, required=True)
    report_parser.add_argument(
        "--llvm-cov",
        help="override the llvm-cov selected in the build's CMake cache",
    )
    report_parser.add_argument(
        "--llvm-profdata",
        help="override the llvm-profdata selected in the build's CMake cache",
    )
    census_parser = subparsers.add_parser("census", help="print the lexical site census")
    census_parser.add_argument("--source-root", type=Path, default=Path.cwd())
    census_parser.add_argument("--build-dir", type=Path)
    subparsers.add_parser("self-test", help="run deterministic lexer/segment tests")
    return parser


def main() -> int:
    args = _parser().parse_args()
    try:
        if args.command == "prepare":
            prepare(args.build_dir, args.source_root)
        elif args.command == "report":
            report(
                args.source_root,
                args.build_dir,
                args.exceptions,
                args.json_output,
                args.markdown_output,
                _coverage_tool(
                    args.build_dir,
                    args.llvm_cov,
                    "SLIDE_LLVM_COV_EXECUTABLE",
                ),
                _coverage_tool(
                    args.build_dir,
                    args.llvm_profdata,
                    "SLIDE_LLVM_PROFDATA_EXECUTABLE",
                ),
            )
        elif args.command == "census":
            census(args.source_root, args.build_dir)
        else:
            self_test()
    except (CoverageError, OSError, json.JSONDecodeError) as error:
        print(f"status-failure coverage: ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
