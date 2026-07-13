"""Report stage of the status-failure coverage gate."""

from __future__ import annotations
import datetime as dt
import json
from pathlib import Path
import tempfile

from .constants import REPORT_SCHEMA
from .scanner import (
    CoverageError,
    Site,
    _matching_statement,
    _status_sequences,
    _tokens,
    scan_sites,
)
from .session import (
    _coverage_tool,
    _load_session,
    _merge_profiles,
    _profile_groups,
    _run_anchor,
    _tool_major,
    read_manifest,
)
from .classify import (
    _counts_for_sites,
    _coverage_files,
    _export,
    _load_exceptions,
    _segment_interval_count,
    classify_active_sites,
)

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
