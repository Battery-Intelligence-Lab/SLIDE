"""Classify stage of the status-failure coverage gate."""

from __future__ import annotations
import bisect
import json
from pathlib import Path
import re
import shlex
import subprocess
from typing import Sequence

from .constants import FAILURE_CLASSES
from .scanner import (
    CoverageError,
    Site,
)
from .session import (
    _run,
)

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
