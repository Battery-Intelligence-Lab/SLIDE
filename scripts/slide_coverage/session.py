"""Session stage of the status-failure coverage gate."""

from __future__ import annotations
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time
from typing import Iterable, Sequence

from .constants import SESSION_SCHEMA, SOURCE_SUFFIXES
from .scanner import (
    CoverageError,
    ManifestEntry,
    _sha256,
)

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
