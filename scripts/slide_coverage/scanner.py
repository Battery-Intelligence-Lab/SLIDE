"""Scanner stage of the status-failure coverage gate."""

from __future__ import annotations
import bisect
import dataclasses
import hashlib
from pathlib import Path
import re
from typing import Sequence
from .constants import SOURCE_SUFFIXES

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
