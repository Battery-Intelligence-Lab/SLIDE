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
import json
from pathlib import Path
import sys


FAILURE_CLASSES = {
    "defensive-only",
    "platform-specific",
    "unreachable-by-construction",
}
SOURCE_SUFFIXES = {".cc", ".cpp", ".cu", ".h", ".hh", ".hpp"}
SESSION_SCHEMA = 1
REPORT_SCHEMA = 1

from slide_coverage.scanner import (
    CoverageError,
    scan_sites,
)
from slide_coverage.session import (
    _coverage_tool,
    prepare,
)
from slide_coverage.classify import (
    classify_active_sites,
)
from slide_coverage.report import (
    report,
    self_test,
)

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
