"""Extract and execute the runnable quickstarts embedded in the v4 docs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote


MARKER = re.compile(r"<!--\s*doc-test:(cpp|python|matlab)\s*-->")
FENCE = re.compile(r"```([A-Za-z0-9_+-]+)\s*$")
EXTENSION = {"cpp": ".cpp", "python": ".py", "matlab": ".m"}
LANGUAGE = {"cpp": "cpp", "python": "python", "matlab": "matlab"}
LINK = re.compile(r"(?<!!)\[[^]]+\]\(([^)\s]+)(?:\s+['\"][^'\"]*['\"])?\)")
HEADING = re.compile(r"^#{1,6}\s+(.+?)\s*$")


def validate_front_matter(markdown: Path) -> None:
    """Require the layout metadata needed for a themed Jekyll page."""

    lines = markdown.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "---":
        raise RuntimeError(f"{markdown}: missing YAML front matter")
    try:
        closing = lines.index("---", 1)
    except ValueError as error:
        raise RuntimeError(f"{markdown}: unterminated YAML front matter") from error
    fields = {
        key.strip(): value.strip()
        for line in lines[1:closing]
        if ":" in line
        for key, value in [line.split(":", 1)]
    }
    if fields.get("layout") != "default":
        raise RuntimeError(f"{markdown}: Jekyll page must use layout: default")
    if not fields.get("title"):
        raise RuntimeError(f"{markdown}: Jekyll page must declare a title")


def heading_ids(markdown: Path) -> set[str]:
    """Return the generated IDs for ATX headings used by the v4 pages."""

    ids: set[str] = set()
    counts: dict[str, int] = {}
    for line in markdown.read_text(encoding="utf-8").splitlines():
        match = HEADING.match(line)
        if match is None:
            continue
        title = re.sub(r"[`*_~]", "", match.group(1)).strip().lower()
        base = re.sub(r"[^\w\s-]", "", title)
        base = re.sub(r"[\s-]+", "-", base).strip("-")
        occurrence = counts.get(base, 0)
        counts[base] = occurrence + 1
        ids.add(base if occurrence == 0 else f"{base}-{occurrence}")
    return ids


def validate_local_links(repo: Path) -> None:
    """Fail when a user-facing v4 Markdown link has no local target or heading."""

    site_sources = [repo / "docs" / "index.md"]
    site_sources.extend(sorted((repo / "docs" / "v4").glob("*.md")))
    for source in site_sources:
        validate_front_matter(source)
    sources = [repo / "README.md", *site_sources]
    checked = 0
    for source in sources:
        text = source.read_text(encoding="utf-8")
        for raw_target in LINK.findall(text):
            target = unquote(raw_target.strip("<>"))
            if target.startswith(("https://", "http://", "mailto:")):
                continue
            path_text, separator, fragment = target.partition("#")
            candidate = source if not path_text else source.parent / path_text
            if candidate.suffix == ".html":
                candidate = candidate.with_suffix(".md")
            candidate = candidate.resolve()
            if not candidate.exists():
                raise RuntimeError(f"{source}: broken local link {raw_target!r}")
            if separator and fragment and candidate.suffix.lower() == ".md":
                if fragment not in heading_ids(candidate):
                    raise RuntimeError(
                        f"{source}: missing heading #{fragment} in {candidate}"
                    )
            checked += 1
    print(f"validated {checked} local documentation links", flush=True)


def extract_blocks(docs: Path) -> dict[str, str]:
    """Return the one marked fenced block for every supported language."""

    blocks: dict[str, str] = {}
    for markdown in sorted(docs.rglob("*.md")):
        lines = markdown.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            marker = MARKER.fullmatch(line.strip())
            if marker is None:
                continue
            name = marker.group(1)
            if name in blocks:
                raise RuntimeError(f"duplicate doc-test:{name} marker")
            if index + 1 >= len(lines):
                raise RuntimeError(f"{markdown}: marker has no fenced block")
            opening = FENCE.fullmatch(lines[index + 1].strip())
            if opening is None or opening.group(1) != LANGUAGE[name]:
                raise RuntimeError(
                    f"{markdown}: doc-test:{name} must precede a {LANGUAGE[name]} fence"
                )
            body: list[str] = []
            for candidate in lines[index + 2 :]:
                if candidate.strip() == "```":
                    blocks[name] = "\n".join(body) + "\n"
                    break
                body.append(candidate)
            else:
                raise RuntimeError(f"{markdown}: unterminated doc-test:{name} fence")

    missing = sorted(set(EXTENSION) - set(blocks))
    if missing:
        raise RuntimeError(f"missing documented quickstarts: {', '.join(missing)}")
    return blocks


def run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    """Run one validation command and fail with its original exit status."""

    print("+", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def test_cpp(repo: Path, output: Path, source: Path) -> None:
    """Build the exact C++ block as an external CMake consumer and run it."""

    project = output / "cpp-project"
    project.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, project / "quickstart.cpp")
    cmake_text = f"""\
cmake_minimum_required(VERSION 3.31)
project(SLIDEDocsQuickstart LANGUAGES CXX)
set(SLIDE_CORE_ONLY ON CACHE BOOL \"\" FORCE)
set(SLIDE_BUILD_CORE_TESTS OFF CACHE BOOL \"\" FORCE)
set(SLIDE_BUILD_PYTHON OFF CACHE BOOL \"\" FORCE)
set(SLIDE_WITH_CUDA OFF CACHE BOOL \"\" FORCE)
set(SLIDE_WITH_MATLAB OFF CACHE BOOL \"\" FORCE)
set(SLIDE_WITH_ZSTD OFF CACHE BOOL \"\" FORCE)
set(SLIDE_WITH_ARROW OFF CACHE BOOL \"\" FORCE)
add_subdirectory(\"{repo.as_posix()}\" slide)
add_executable(slide_docs_quickstart quickstart.cpp)
target_link_libraries(slide_docs_quickstart PRIVATE slide_core)
enable_testing()
add_test(NAME slide_docs_quickstart COMMAND slide_docs_quickstart)
"""
    (project / "CMakeLists.txt").write_text(cmake_text, encoding="utf-8")
    build = output / "cpp-build"
    run(
        [
            "cmake",
            "-S",
            str(project),
            "-B",
            str(build),
            "-DCMAKE_BUILD_TYPE=Release",
        ],
        cwd=repo,
    )
    run(
        ["cmake", "--build", str(build), "--config", "Release", "--parallel", "2"],
        cwd=repo,
    )
    run(
        [
            "ctest",
            "--test-dir",
            str(build),
            "-C",
            "Release",
            "--output-on-failure",
        ],
        cwd=repo,
    )


def test_python(repo: Path, source: Path) -> None:
    """Run the exact Python block away from the source tree."""

    run([sys.executable, str(source)], cwd=source.parent)


def test_matlab(repo: Path, source: Path, executable: str | None) -> None:
    """Run the exact MATLAB block with the repository package on its path."""

    matlab = executable or shutil.which("matlab")
    if matlab is None:
        raise RuntimeError("MATLAB requested but no executable was supplied or found")
    env = os.environ.copy()
    env["SLIDE_ROOT"] = str(repo)
    escaped = source.resolve().as_posix().replace("'", "''")
    run([matlab, "-batch", f"run('{escaped}')"], cwd=repo, env=env)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--output", type=Path, required=True)
    execution = parser.add_mutually_exclusive_group(required=True)
    execution.add_argument(
        "--language",
        action="append",
        choices=tuple(EXTENSION),
        help="quickstart to execute; repeat for multiple languages",
    )
    execution.add_argument(
        "--extract-only",
        action="store_true",
        help="validate and extract every quickstart without executing one",
    )
    parser.add_argument("--matlab", help="path to the MATLAB executable")
    args = parser.parse_args()

    repo = args.repo.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    validate_local_links(repo)
    blocks = extract_blocks(repo / "docs")
    sources: dict[str, Path] = {}
    for name, text in blocks.items():
        source = output / f"quickstart{EXTENSION[name]}"
        source.write_text(text, encoding="utf-8")
        sources[name] = source

    for language in dict.fromkeys(args.language or []):
        if language == "cpp":
            test_cpp(repo, output, sources[language])
        elif language == "python":
            test_python(repo, sources[language])
        else:
            test_matlab(repo, sources[language], args.matlab)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
