"""Make fortran case consistent and reformat for more consistent spacing.

Run without arguments, or with ``--all``, to update every tracked free-form
Fortran source in this repository (``.f90``, ``.f95``, ``.f03``, ``.f08``,
``.f18``, or ``.f23``, including uppercase variants). Positional paths limit
the update to those files.

Quoted text and C-preprocessor lines are left unchanged; comment text is preserved
apart from normalizing the comment marker's surrounding spaces.
Code lines longer than 120 columns are re-wrapped by default using free-form
Fortran continuations; use --no-wrap to disable re-wrapping.
Macros supplied by the build can be declared with ``-D NAME[=VALUE]`` or
``--define NAME[=VALUE]``; that spelling becomes the formatter's canonical macro case.
Commented-out assignment statements have their assignment and addition operators normalized.
OpenMP directive keywords are normalized to uppercase.
"""

import argparse
import difflib
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


_GIT_HOOK_CONTEXT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR", "GIT_INDEX_FILE")


def _git_env() -> dict[str, str]:
    """Environment for nested ``git`` calls with any enclosing hook's repository context removed.

    When this script runs as a git hook, git sets ``GIT_DIR`` (but not ``GIT_WORK_TREE``)
    in the environment. A nested ``git`` invocation that inherits it reports its own
    ``cwd`` as the toplevel instead of searching upward for the real one, so these must
    not be passed through.
    """
    return {key: value for key, value in os.environ.items() if key not in _GIT_HOOK_CONTEXT_VARS}


def find_repository_root() -> Path | None:
    """Return the enclosing git checkout for the script, if there is one."""
    try:
        result = subprocess.run(
            ["git", "-C", str(SCRIPT_DIR), "rev-parse", "--show-toplevel"],
            env=_git_env(),
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return None
    if result.returncode == 0 and result.stdout.strip():
        return Path(result.stdout.strip()).resolve()
    return None


REPOSITORY_ROOT = find_repository_root()
IDENTIFIER = r"[A-Za-z][A-Za-z0-9_]*"
USE_MODULE = re.compile(rf"^\s*use\b(?:\s*,[^:\n]*::\s*|\s*::\s*|\s+)({IDENTIFIER})\b", re.IGNORECASE)
TOKEN = re.compile(r"\b[A-Za-z][A-Za-z0-9_]*\b")
REAL_LITERAL_EXPONENT = re.compile(
    r"(?<![A-Za-z0-9_])(?P<mantissa>(?:\d+(?:\.\d*)?|\.\d+))(?P<marker>[EeDd])(?P<exponent>[+-]?\d+)"
)
DIRECTIVE_SENTINEL = re.compile(r"!(?:\$|(?:dir|dec|gcc)\$)", re.IGNORECASE)
LEGACY_OPERATOR = re.compile(r"\.(eq|ne|lt|le|gt|ge)\.", re.IGNORECASE)
MODERN_OPERATOR = {
    "eq": "==",
    "ne": "/=",
    "lt": "<",
    "le": "<=",
    "gt": ">",
    "ge": ">=",
}
# A bare `<` or `>` next to another angle bracket, or right after `-`, is not a
# Fortran operator: it is prose in a comment (`k*tau>>1`, `beta->1`), and
# spacing it out would rewrite the text.
SPACED_OPERATOR = re.compile(
    r"=>|==|/=|<=|>=|(?<![=<>])<(?![<>])|(?<![=<>-])>(?![<>])|(?<![<>=/])=(?!=|>)|\.(?:and|or|not|eqv|neqv)\.",
    re.IGNORECASE,
)
ARITHMETIC_OPERATOR = re.compile(r"\*\*|//|[+*/-]")
COMPACT_ARITHMETIC_OPERATORS = frozenset(("*", "/", "**"))
MAX_LINE_LENGTH = 120
FORTRAN_SOURCE_EXTENSIONS = frozenset({".f90", ".f95", ".f03", ".f08", ".f18", ".f23"})
FORTRAN_SOURCE_PATTERNS = tuple(
    f"*{extension}"
    for extension in sorted(FORTRAN_SOURCE_EXTENSIONS | {ext.upper() for ext in FORTRAN_SOURCE_EXTENSIONS})
)
COMPOUND_KEYWORDS = {
    "blockdata": "block data",
    "elseif": "else if",
    "endif": "end if",
    "endassociate": "end associate",
    "endblock": "end block",
    "endblockdata": "end block data",
    "enddo": "end do",
    "endfile": "end file",
    "endforall": "end forall",
    "endfunction": "end function",
    "endinterface": "end interface",
    "endmodule": "end module",
    "endprogram": "end program",
    "endselect": "end select",
    "endsubroutine": "end subroutine",
    "endtype": "end type",
    "endwhere": "end where",
}
GO_TO = re.compile(r"(?<![A-Za-z0-9_])go[ \t]+to\b", re.IGNORECASE)
MULTIWORD_KEYWORD_PAIRS = (
    ("abstract", "interface"),
    ("change", "team"),
    ("class", "default"),
    ("class", "is"),
    ("do", "concurrent"),
    ("double", "precision"),
    ("event", "post"),
    ("event", "wait"),
    ("fail", "image"),
    ("form", "team"),
    ("impure", "elemental"),
    ("pure", "elemental"),
    ("rank", "default"),
    ("select", "case"),
    ("select", "rank"),
    ("sync", "all"),
    ("sync", "images"),
    ("sync", "memory"),
    ("sync", "team"),
    ("type", "default"),
    ("type", "is"),
)
MULTIWORD_KEYWORD_SPACING = re.compile(
    r"\b(?:" + "|".join(rf"{first}[ \t]+{second}" for first, second in MULTIWORD_KEYWORD_PAIRS) + r")\b",
    re.IGNORECASE,
)

COMPOUND_KEYWORD = re.compile(
    rf"^([ \t]*)(?P<keyword>{'|'.join(sorted(COMPOUND_KEYWORDS, key=len, reverse=True))})\b",
    re.IGNORECASE,
)
COMMON_BLOCK_PREFIX = re.compile(rf"^([ \t]*common)\s*/\s*({IDENTIFIER})\s*/\s*", re.IGNORECASE)
PARENTHESIZED_STATEMENT_NAMES = frozenset(
    {
        "allocate",
        "allocated",
        "backspace",
        "close",
        "deallocate",
        "endfile",
        "flush",
        "inquire",
        "nullify",
        "open",
        "read",
        "rewind",
        "wait",
        "write",
    }
)
PARENTHESIZED_STATEMENT = re.compile(
    r"(?<![A-Za-z0-9_])"
    rf"(?:{'|'.join(sorted(PARENTHESIZED_STATEMENT_NAMES))})"
    r"[ \t]*\(",
    re.IGNORECASE,
)
EMPTY_SUBROUTINE_ARGUMENTS = re.compile(
    rf"^(\s*.*?\bsubroutine\s+{IDENTIFIER})\s*\(\s*\)",
    re.IGNORECASE,
)
DECLARATION_ATTRIBUTES = frozenset(
    {
        "allocatable",
        "asynchronous",
        "codimension",
        "contiguous",
        "dimension",
        "external",
        "intent",
        "intrinsic",
        "optional",
        "parameter",
        "pointer",
        "private",
        "protected",
        "public",
        "save",
        "target",
        "value",
        "volatile",
    }
)

# Break-point candidates for line wrapping. wrap_position() prefers the
# shallowest bracket-nesting depth first, then the earliest tier below, so a
# wrap lands on the same outermost boundary a human would pick instead of
# splitting mid-expression. Every pattern matches a whole Fortran token, since
# a break placed inside one (`/` of `//`, `/` of `/=`, `/` of `/)`) would need
# a leading `&` on the continuation line to stay valid.
# Assignment (`=`/`=>`) is deliberately excluded from this tiered search: it
# sits at depth 0 right after the LHS, so ranking it by depth would make it
# beat RHS operators nested inside any parentheses (e.g. `x = f(a + b)`),
# even though those inner operators usually pack the line far better. It is
# handled separately in wrap_position(), as the end of the statement head.
COMMA_BREAK_OPERATOR = re.compile(r",")
EQUIVALENCE_BREAK_OPERATOR = re.compile(r"\.(?:eqv|neqv)\.", re.IGNORECASE)
DISJUNCTION_BREAK_OPERATOR = re.compile(r"\.or\.", re.IGNORECASE)
CONJUNCTION_BREAK_OPERATOR = re.compile(r"\.and\.", re.IGNORECASE)
COMPARISON_BREAK_OPERATOR = re.compile(r"==|/=|<=|>=|(?<![=])<|(?<![=])>")
CONCATENATION_BREAK_OPERATOR = re.compile(r"//")
ADDITIVE_BREAK_OPERATOR = re.compile(r"(?<=\s)[+-](?=\s)")
MULTIPLICATIVE_BREAK_OPERATOR = re.compile(r"\*\*|\*|(?<!/)/(?![=)/])")
# Ordered from loosest to tightest binding: a comma separates whole items, and
# each following tier binds its operands more tightly than the one before it.
BREAK_OPERATORS: tuple[re.Pattern[str], ...] = (
    COMMA_BREAK_OPERATOR,
    EQUIVALENCE_BREAK_OPERATOR,
    DISJUNCTION_BREAK_OPERATOR,
    CONJUNCTION_BREAK_OPERATOR,
    COMPARISON_BREAK_OPERATOR,
    CONCATENATION_BREAK_OPERATOR,
    ADDITIVE_BREAK_OPERATOR,
    MULTIPLICATIVE_BREAK_OPERATOR,
)
ASSIGNMENT_BREAK_OPERATOR = re.compile(r"=>|(?<![<>=/])=(?!=|>)")
# Only treat a comment as commented-out code when its first non-whitespace
# text is an assignable-looking name. An equals sign elsewhere is prose.
COMMENTED_ASSIGNMENT = re.compile(
    rf"^![ \t]*{IDENTIFIER}(?:[ \t]*(?:%[ \t]*{IDENTIFIER}|\([^!\n]*\)))*[ \t]*(?<![<>=/])=(?!=|>)"
)
# Share of the available width a break must use before it counts as filling its
# line. An outermost operator can sit two characters in (`d*(...)`), where
# taking it would strand one operand on a line of its own; such a break is
# ranked behind every well-filled one and used only if nothing else is left.
MINIMUM_BREAK_FILL = 0.25

# Fortran 90 through 2018 language keywords, including the multi-word
# statement components that are normalized independently. Keep identifiers
# out of this list even if they happen to be uppercase in the current source.
FORTRAN_KEYWORDS = frozenset(
    [
        "all",
        "abstract",
        "allocatable",
        "allocate",
        "assign",
        "assigned",
        "assignment",
        "associate",
        "asynchronous",
        "backspace",
        "bind",
        "block",
        "blockdata",
        "call",
        "case",
        "character",
        "class",
        "change",
        "close",
        "codimension",
        "common",
        "complex",
        "concurrent",
        "contains",
        "continue",
        "contiguous",
        "critical",
        "cycle",
        "data",
        "deallocate",
        "default",
        "deferred",
        "dimension",
        "do",
        "double",
        "else",
        "elseif",
        "elsewhere",
        "elemental",
        "end",
        "endassociate",
        "endblock",
        "endblockdata",
        "enddo",
        "endfile",
        "endforall",
        "endfunction",
        "endif",
        "endinterface",
        "endmodule",
        "endprogram",
        "endselect",
        "endsubroutine",
        "endtype",
        "endwhere",
        "entry",
        "enum",
        "enumerator",
        "equivalence",
        "error",
        "exit",
        "extends",
        "external",
        "event",
        "fail",
        "final",
        "flush",
        "forall",
        "format",
        "form",
        "function",
        "generic",
        "go",
        "to",
        "goto",
        "if",
        "implicit",
        "import",
        "in",
        "include",
        "inquire",
        "integer",
        "intent",
        "interface",
        "intrinsic",
        "kind",
        "image",
        "images",
        "impure",
        "logical",
        "lock",
        "local",
        "local_init",
        "memory",
        "module",
        "mold",
        "namelist",
        "none",
        "non_intrinsic",
        "non_overridable",
        "non_recursive",
        "nopass",
        "nullify",
        "only",
        "open",
        "operator",
        "post",
        "optional",
        "out",
        "inout",
        "parameter",
        "pass",
        "pause",
        "pointer",
        "print",
        "private",
        "procedure",
        "program",
        "protected",
        "public",
        "pure",
        "rank",
        "read",
        "real",
        "recursive",
        "reduce",
        "result",
        "return",
        "rewrite",
        "rewind",
        "save",
        "select",
        "sequence",
        "shared",
        "source",
        "stop",
        "submodule",
        "subroutine",
        "sync",
        "team",
        "target",
        "then",
        "type",
        "unlock",
        "until_count",
        "use",
        "value",
        "volatile",
        "wait",
        "where",
        "while",
        "write",
        "notify",
    ]
)

# Standard statement specifiers are not reserved words, but consistently
# lowercasing them in their ordinary spelling makes keyword arguments such as
# ``source=`` and ``mold=`` match the rest of the formatted language. Local
# declarations still take precedence in ``lowercase_keyword``.
FORTRAN_SPECIFIERS = frozenset(
    {
        "acquired_lock",
        "access",
        "action",
        "advance",
        "blank",
        "decimal",
        "delim",
        "direct",
        "encoding",
        "eor",
        "err",
        "errmsg",
        "exist",
        "file",
        "fmt",
        "form",
        "formatted",
        "id",
        "iomsg",
        "iostat",
        "newunit",
        "nextrec",
        "nml",
        "number",
        "name",
        "new_index",
        "opened",
        "pad",
        "pending",
        "pos",
        "position",
        "quiet",
        "readwrite",
        "rec",
        "recl",
        "round",
        "sequential",
        "sign",
        "size",
        "stat",
        "status",
        "stream",
        "leading_zero",
        "unformatted",
        "unit",
    }
)
FORTRAN_STANDARD_WORDS = FORTRAN_KEYWORDS | FORTRAN_SPECIFIERS

# Intrinsic procedures are globally available lower-case names, rather than
# language keywords. A declaration in the current scope can therefore shadow
# one without being rewritten to lower case.
INTRINSIC_PROCEDURES = frozenset(
    {
        "abs",
        "acos",
        "allocated",
        "asin",
        "atan",
        "atan2",
        "ceiling",
        "cmplx",
        "conjg",
        "cos",
        "cosh",
        "cpu_time",
        "dim",
        "dot_product",
        "exp",
        "floor",
        "huge",
        "iand",
        "ibclr",
        "ibits",
        "ibset",
        "ieor",
        "index",
        "int",
        "ishft",
        "iso_fortran_env",
        "is_iostat_end",
        "is_iostat_eor",
        "lbound",
        "len",
        "len_trim",
        "log",
        "log10",
        "max",
        "maxloc",
        "maxval",
        "merge",
        "min",
        "minloc",
        "minval",
        "mod",
        "modulo",
        "nint",
        "pack",
        "precision",
        "product",
        "random_number",
        "repeat",
        "reshape",
        "sign",
        "sin",
        "sinh",
        "size",
        "sqrt",
        "tan",
        "tanh",
        "tiny",
        "trim",
        "ubound",
        "unpack",
        "verify",
    }
)

# OpenMP directive and clause vocabulary is uppercased only inside OpenMP
# directives; these names are not Fortran keywords in ordinary code.
OPENMP_KEYWORDS = frozenset(
    {
        "omp",
        "do",
        "atomic",
        "barrier",
        "cancel",
        "cancellation",
        "critical",
        "declare",
        "distribute",
        "end",
        "flush",
        "loop",
        "master",
        "masked",
        "ordered",
        "parallel",
        "sections",
        "section",
        "simd",
        "single",
        "target",
        "task",
        "taskgroup",
        "taskloop",
        "taskwait",
        "taskyield",
        "teams",
        "threadprivate",
        "workshare",
        "allocate",
        "collapse",
        "copyin",
        "copyprivate",
        "default",
        "firstprivate",
        "if",
        "lastprivate",
        "linear",
        "map",
        "nowait",
        "num_threads",
        "private",
        "reduction",
        "schedule",
        "static",
        "dynamic",
        "guided",
        "runtime",
        "shared",
        "simdlen",
        "proc_bind",
        "defaultmap",
        "depend",
        "device",
        "dist_schedule",
        "final",
        "grainsize",
        "hint",
        "in_reduction",
        "is_device_ptr",
        "mergeable",
        "nogroup",
        "num_tasks",
        "order",
        "priority",
        "safelen",
        "thread_limit",
        "to",
        "from",
        "use_device_addr",
        "use_device_ptr",
    }
)

INTRINSIC_NAMES = INTRINSIC_PROCEDURES | frozenset(
    {"false", "true", "and", "or", "not", "eq", "ne", "lt", "le", "gt", "ge", "eqv", "neqv"}
)
CPP_DEFINE = re.compile(rf"^\s*#\s*define\s+({IDENTIFIER})\b")


def is_preprocessor_line(line: str) -> bool:
    return line.lstrip().startswith("#")


def cpp_line_continues(line: str) -> bool:
    return line.rstrip("\r\n").endswith("\\")


def extract_preprocessor_cases(source: str) -> dict[str, str]:
    """Return unambiguous case-sensitive C-preprocessor macro spellings."""
    occurrences: dict[str, set[str]] = {}
    for line in source.splitlines():
        match = CPP_DEFINE.match(line)
        if match:
            name = match.group(1)
            occurrences.setdefault(name.lower(), set()).add(name)
    return {key: next(iter(values)) for key, values in occurrences.items() if len(values) == 1}


def extract_preprocessor_names(source: str) -> frozenset[str]:
    """Return case-sensitive C-preprocessor macro names used by this source."""
    return frozenset(extract_preprocessor_cases(source).values())


def normalize_macro_cases(macros: Mapping[str, str] | Collection[str] = ()) -> dict[str, str]:
    """Normalize explicit macro names to a case-insensitive lookup mapping."""
    if isinstance(macros, Mapping):
        return {str(key).lower(): str(value) for key, value in macros.items()}
    return {name.lower(): name for name in macros}


def replace_preprocessor_cases(source: str, macro_cases: Mapping[str, str]) -> str:
    """Apply requested macro spelling to unquoted Fortran tokens, never CPP bodies."""
    if not macro_cases:
        return source
    output: list[str] = []
    _, states = _scan_source(source)
    for state in states:
        line = state.text
        if state.is_cpp:
            output.append(line)
            continue
        line_output: list[str] = []
        quote = state.quote_in
        index = 0
        while index < len(line):
            char = line[index]
            if quote:
                line_output.append(char)
                if char == quote:
                    if line[index + 1 : index + 2] == quote:
                        line_output.append(line[index + 1])
                        index += 2
                        continue
                    quote = None
                index += 1
                continue
            if char in "\"'":
                quote = char
                line_output.append(char)
                index += 1
                continue
            if char == "!":
                line_output.append(line[index:])
                break
            match = TOKEN.match(line, index)
            if match:
                token = match.group()
                line_output.append(macro_cases.get(token.lower(), token))
                index = match.end()
            else:
                line_output.append(char)
                index += 1
        output.append("".join(line_output))
    return "".join(output)


@dataclass(frozen=True)
class DeclaredName:
    """A case-sensitive spelling found in a Fortran declaration."""

    kind: str
    name: str


@dataclass(frozen=True)
class FileDeclarationCases:
    """The declaration spellings that are unambiguous for one source file."""

    module_cases: Mapping[str, str]
    symbol_cases: Mapping[str, str]
    procedure_cases: tuple["ProcedureDeclarationCases", ...] = ()
    scope_cases: tuple["NamedScopeCase", ...] = ()
    type_procedure_cases: Mapping[str, str] = field(default_factory=dict)
    type_component_cases: Mapping[tuple[str, str], str] = field(default_factory=dict)
    variable_type_cases: Mapping[str, str] = field(default_factory=dict)
    type_component_type_cases: Mapping[tuple[str, str], str] = field(default_factory=dict)


@dataclass(frozen=True)
class ProcedureDeclarationCases:
    """The local variable spellings for one function or subroutine scope."""

    start_line: int
    end_line: int
    local_cases: Mapping[str, str]
    local_names: frozenset[str] = frozenset()
    local_types: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class NamedScopeCase:
    """The declaration spelling for a named module, program, or procedure scope."""

    kind: str
    start_line: int
    end_line: int
    name: str


@dataclass(frozen=True)
class ScopedDeclaredNames:
    """Type, procedure, and module-variable names declared within one module,
    program, function, or subroutine's textual body.

    Grouping by enclosing scope, rather than collecting one flat set for the
    whole file, keeps a declaration made in one module -- public or private
    -- from leaking into an unrelated module later in the same file: a name
    is only trustworthy for overriding keyword/intrinsic lowercasing where it
    is actually in scope.
    """

    start_line: int
    end_line: int
    names: frozenset[str]


@dataclass(frozen=True)
class CodeStatement:
    """A code statement and its inclusive physical source line range."""

    start_line: int
    end_line: int
    text: str


@dataclass(frozen=True)
class PhysicalLine:
    """Lexical state for one physical source line.

    Centralizing this state prevents formatter passes from disagreeing about
    preprocessor continuations, strings, comments, and Fortran continuations.
    """

    number: int
    text: str
    is_cpp: bool
    cpp_continuation_in: bool
    cpp_continuation_out: bool
    is_blank: bool
    is_comment: bool
    quote_in: str | None
    quote_out: str | None
    continuation_in: bool
    continuation_out: bool


def scan_physical_lines(lines: Iterable[str]) -> tuple[PhysicalLine, ...]:
    """Return shared lexical/continuation state for *lines*."""
    states: list[PhysicalLine] = []
    cpp_continuation = False
    quote: str | None = None
    continuation = False
    for number, line in enumerate(lines):
        cpp_in = cpp_continuation
        is_cpp = cpp_in or is_preprocessor_line(line)
        is_blank = not line.strip()
        is_comment = not is_cpp and line.lstrip().startswith("!")
        quote_in = quote
        continuation_in = continuation
        if is_cpp:
            cpp_continuation = cpp_line_continues(line)
            quote_out = None
            continuation_out = False
            quote = None
            continuation = False
        elif is_blank or is_comment:
            cpp_continuation = False
            # Comments/blank lines may legally appear between continued source
            # lines. They do not terminate the pending statement or string.
            quote_out = quote
            continuation_out = continuation
        else:
            cpp_continuation = False
            continuation_out, quote_out = continues_statement(line, quote)
            quote = quote_out
            continuation = continuation_out
        states.append(
            PhysicalLine(
                number,
                line,
                is_cpp,
                cpp_in,
                cpp_continuation,
                is_blank,
                is_comment,
                quote_in,
                quote_out,
                continuation_in,
                continuation_out,
            )
        )
    return tuple(states)


@lru_cache(maxsize=128)
def _scan_source(source: str) -> tuple[tuple[str, ...], tuple[PhysicalLine, ...]]:
    """Split and lex *source* once, reusing physical-line state across passes."""
    lines = tuple(source.splitlines(keepends=True))
    return lines, scan_physical_lines(lines)


def _statement_end_index(states: tuple[PhysicalLine, ...], start: int) -> int:
    """Return the last physical line belonging to a Fortran statement."""
    state = states[start]
    if state.is_cpp or state.is_blank or state.is_comment or not state.continuation_out:
        return start
    continuation = True
    index = start
    while continuation and index + 1 < len(states):
        index += 1
        current = states[index]
        if current.is_cpp or current.is_blank or current.is_comment:
            continue
        continuation = current.continuation_out
    return index


def _join_statement_parts(parts: Iterable[tuple[str, bool]]) -> str:
    """Join physical pieces, omitting a space only for split-token continuations."""
    joined = ""
    for part, joins_previous in parts:
        if not part:
            continue
        if joined and joins_previous:
            joined += part
        elif joined:
            joined += " " + part
        else:
            joined = part
    return joined


def _statement_parts(statement: Iterable[str]) -> list[tuple[str, bool]]:
    """Remove continuation markers while retaining whether each boundary joins tokens."""
    lines = list(statement)
    contexts = [code_context(line) for line in lines]
    parts: list[tuple[str, bool]] = []
    for index, context in enumerate(contexts):
        code = context.strip()
        leading = context.lstrip(" \t")
        leading_token_continuation = (
            index > 0 and leading.startswith("&") and len(leading) > 1 and not leading[1].isspace()
        )
        if index:
            code = re.sub(r"^&\s*", "", code)
        if index + 1 < len(lines):
            code = re.sub(r"\s*&$", "", code)
        if code:
            previous = contexts[index - 1].rstrip("\r\n \t") if index else ""
            joins_previous = (
                leading_token_continuation
                and previous.endswith("&")
                and len(previous) > 1
                and not previous[-2].isspace()
            )
            parts.append((code, joins_previous))
    return parts


def _trailing_token_continuation(line: str) -> bool:
    """Return whether a continuation marker touches the preceding token."""
    code = code_context(line).rstrip("\r\n \t")
    return code.endswith("&") and len(code) > 1 and not code[-2].isspace()


@lru_cache(maxsize=128)
def _code_statements(source: str) -> tuple[CodeStatement, ...]:
    """Return logical Fortran statements, reusing the result across analysis passes."""
    lines, states = _scan_source(source)
    statements: list[CodeStatement] = []
    index = 0
    while index < len(lines):
        state = states[index]
        if state.is_cpp or state.is_blank or state.is_comment:
            index += 1
            continue
        end_index = _statement_end_index(states, index)
        if end_index == index:
            text = code_context(lines[index]).strip()
        else:
            statement_lines = (
                lines[pos]
                for pos in range(index, end_index + 1)
                if not states[pos].is_cpp and not states[pos].is_blank and not states[pos].is_comment
            )
            parts = _statement_parts(statement_lines)
            text = _join_statement_parts(parts) if parts else ""
        if text:
            for statement_text in _split_top_level_statements(text):
                if statement_text:
                    statements.append(CodeStatement(index, end_index, statement_text))
        index = end_index + 1
    return tuple(statements)


def _iter_code_statements_with_lines(source: str) -> Iterable[CodeStatement]:
    """Yield Fortran statements while ignoring CPP and interleaved comments."""
    return iter(_code_statements(source))


def _iter_code_statements(source: str) -> Iterable[str]:
    """Yield continued Fortran statements with comments and layout removed."""
    return (statement.text for statement in _code_statements(source))


DECLARED_TYPE_HEADER = re.compile(r"^\s*(?:abstract\s+)?type\b(?P<rest>.*)$", re.IGNORECASE)
TYPE_CLASS_CONTEXT = re.compile(r"^\s*(?:type|class)\s+(?:is|default|\()", re.IGNORECASE)
PROCEDURE_WORD = re.compile(r"\b(?:function|subroutine)\b", re.IGNORECASE)
DECLARATION_ENTITY = re.compile(rf"\s*({IDENTIFIER})\b")
DECLARED_VARIABLE_TYPE = re.compile(rf"^\s*(?:type|class)\s*\(\s*({IDENTIFIER})\s*\)", re.IGNORECASE)


def _declared_type_name(statement: str) -> str | None:
    """Return a derived-type name from a type declaration, if present."""
    match = DECLARED_TYPE_HEADER.match(statement)
    if match is None:
        return None
    rest = match.group("rest").lstrip()
    if rest.startswith("(") or re.match(r"(?:is|default)\b", rest, re.IGNORECASE):
        return None
    if rest.startswith(","):
        separator = rest.find("::")
        if separator < 0:
            return None
        rest = rest[separator + 2 :].lstrip()
    elif rest.startswith("::"):
        rest = rest[2:].lstrip()
    match = DECLARATION_ENTITY.match(rest)
    return None if match is None else match.group(1)


DECLARATION_STATEMENT = re.compile(
    r"^\s*(?:integer|real|double\s+precision|complex|logical|character|type|class|procedure|"
    r"dimension|allocatable|pointer|target|optional|parameter|save|value|volatile|asynchronous|"
    r"contiguous|codimension)\b",
    re.IGNORECASE,
)
OLD_STYLE_DECLARATION = re.compile(
    r"^\s*(?:integer|real|double\s+precision|complex|logical|character|type(?:\s*\([^)]*\))?|"
    r"class(?:\s*\([^)]*\))?)"
    r"(?:\s*\([^)]*\)|\s*\*\s*[A-Za-z0-9]+)?\s+(.+)$",
    re.IGNORECASE,
)


def _split_top_level(text: str, separator: str = ",") -> list[str]:
    """Split *text* at separators outside delimiters and character literals."""
    if separator not in text:
        return [text]
    parts: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    for index, char in enumerate(text):
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char in "([":
            depth += 1
        elif char in ")]":
            depth = max(0, depth - 1)
        elif char == separator and depth == 0:
            parts.append(text[start:index])
            start = index + 1
    parts.append(text[start:])
    return parts


def _split_top_level_statements(text: str) -> list[str]:
    """Split semicolon-separated statements outside parentheses and literals."""
    if ";" not in text:
        return [text.strip()]
    parts: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    index = 0
    while index < len(text):
        char = text[index]
        if quote:
            if char == quote:
                if text[index + 1 : index + 2] == quote:
                    index += 2
                    continue
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == ";" and depth == 0:
            parts.append(text[start:index].strip())
            start = index + 1
        index += 1
    parts.append(text[start:].strip())
    return parts


def _declared_variable_names(statement: str) -> list[str]:
    """Return variable names from a type or attribute declaration statement."""
    if not DECLARATION_STATEMENT.match(statement):
        return []
    if _declared_type_name(statement) or TYPE_CLASS_CONTEXT.match(statement):
        return []
    if PROCEDURE_WORD.search(statement):
        return []

    separator = statement.find("::")
    if separator >= 0:
        entities = statement[separator + 2 :]
    else:
        old_style = OLD_STYLE_DECLARATION.match(statement)
        if old_style is None:
            return []
        entities = old_style.group(1)

    names: list[str] = []
    for entity in _split_top_level(entities):
        match = DECLARATION_ENTITY.match(entity)
        if match:
            names.append(match.group(1))
    return names


def _declared_variable_types(statement: str) -> dict[str, str]:
    """Return names and declared type names for a TYPE/CLASS declaration."""
    match = DECLARED_VARIABLE_TYPE.match(statement)
    if match is None:
        return {}
    return {name.lower(): match.group(1) for name in _declared_variable_names(statement)}


def extract_variable_types(source: str) -> dict[str, str]:
    """Return unambiguous declared variable types found in *source*."""
    occurrences: dict[str, set[str]] = {}
    for statement in _iter_code_statements(source):
        for name, type_name in _declared_variable_types(statement).items():
            occurrences.setdefault(name, set()).add(type_name)
    return {name: next(iter(type_names)) for name, type_names in occurrences.items() if len(type_names) == 1}


SCOPE_END_START = re.compile(r"^\s*end\b", re.IGNORECASE)
SCOPE_HEADER = re.compile(rf"\b(function|subroutine|program)\s+({IDENTIFIER})\b", re.IGNORECASE)
RESULT_CLAUSE = re.compile(rf"\bresult\s*\(\s*({IDENTIFIER})\s*\)", re.IGNORECASE)
LEADING_IDENTIFIER = re.compile(rf"\s*({IDENTIFIER})\b")
PROCEDURE_UNIT_END = re.compile(r"^\s*end\s*(?:function|subroutine|program)\b", re.IGNORECASE)
CONTAINS_STATEMENT = re.compile(r"^\s*contains\b", re.IGNORECASE)
TYPE_DECLARATION_START = re.compile(r"^\s*(?:abstract\s+)?type\b", re.IGNORECASE)
SELECT_TYPE_ALIAS = re.compile(rf"^\s*select\s+type\s*\(\s*({IDENTIFIER})\s*=>", re.IGNORECASE)
FORMAT_STATEMENT_START = re.compile(r"^\s*(?:\d+\s+)?format\s*\(", re.IGNORECASE)
IF_STATEMENT_START = re.compile(r"^(\s*(?:else\s+)?if\s*\()", re.IGNORECASE)
FUNCTION_OR_SUBROUTINE_WORD = re.compile(r"\b(function|subroutine)\b", re.IGNORECASE)


def _scope_header(statement: str) -> tuple[str, str, list[str]] | None:
    """Return kind, name, and header argument/result names for a program unit."""
    if SCOPE_END_START.match(statement):
        return None
    match = SCOPE_HEADER.search(statement)
    if match is None:
        return None

    kind = match.group(1).lower()
    name = match.group(2)
    names: list[str] = []
    position = match.end()
    while position < len(statement) and statement[position].isspace():
        position += 1
    if kind != "program" and position < len(statement) and statement[position] == "(":
        depth = 0
        close = None
        for index in range(position, len(statement)):
            if statement[index] == "(":
                depth += 1
            elif statement[index] == ")":
                depth -= 1
                if depth == 0:
                    close = index
                    break
        if close is not None:
            for argument in _split_top_level(statement[position + 1 : close]):
                argument_match = LEADING_IDENTIFIER.match(argument)
                if argument_match:
                    names.append(argument_match.group(1))
            position = close + 1

    result = RESULT_CLAUSE.search(statement, position)
    if result:
        names.append(result.group(1))
    return kind, name, names


def _scope_header_names(statement: str) -> tuple[str, list[str]] | None:
    """Return a procedure/program name and header argument/result names, if declared."""
    header = _scope_header(statement)
    return None if header is None else (header[1], header[2])


def _select_type_alias(statement: str) -> str | None:
    """Return the construct-local associate name from a SELECT TYPE alias."""
    match = SELECT_TYPE_ALIAS.match(statement)
    return None if match is None else match.group(1)


NAMED_SCOPE_END = re.compile(
    rf"^\s*end\s*(module|program|function|subroutine)\b(?:\s+({IDENTIFIER}))?",
    re.IGNORECASE,
)
BARE_PROGRAM_UNIT_END = re.compile(r"^\s*end\s*$", re.IGNORECASE)
MODULE_DECLARATION = re.compile(
    rf"^\s*module\s+(?!(?:procedure|subroutine|function)\b)({IDENTIFIER})\b",
    re.IGNORECASE,
)
PROGRAM_UNIT_END = re.compile(r"^\s*end\s+(?:module|program|function|subroutine)\b", re.IGNORECASE)
TYPE_DEFINITION_START = re.compile(rf"^\s*type(?!\s+is\b)(?:\s*,[^:]*)?\s*(?:::)?\s*{IDENTIFIER}\b", re.IGNORECASE)
TYPE_DEFINITION_END = re.compile(r"^\s*end\s*type\b", re.IGNORECASE)
INTERFACE_START = re.compile(r"^\s*(?:abstract\s+)?interface\b", re.IGNORECASE)
INTERFACE_END = re.compile(r"^\s*end\s*interface\b", re.IGNORECASE)
DECLARATION_PROCEDURE_START = re.compile(r"^\s*procedure\b", re.IGNORECASE)
END_STATEMENT = re.compile(r"^\s*end\b", re.IGNORECASE)
END_MODULE = re.compile(r"^\s*end\s*module\b", re.IGNORECASE)
CONTAINS_STATEMENT = re.compile(r"^\s*contains\b", re.IGNORECASE)
PROCEDURE_DECLARATION_NAME = re.compile(rf"\b(function|subroutine)\s+({IDENTIFIER})\b", re.IGNORECASE)


def active_procedure_at(
    procedure_cases: Iterable[ProcedureDeclarationCases], line_number: int
) -> ProcedureDeclarationCases | None:
    """Return the innermost procedure scope containing *line_number*."""
    return next(
        (
            procedure
            for procedure in reversed(tuple(procedure_cases))
            if procedure.start_line <= line_number <= procedure.end_line
        ),
        None,
    )


def _active_procedures_by_line(
    procedure_cases: Iterable[ProcedureDeclarationCases], line_count: int
) -> tuple[ProcedureDeclarationCases | None, ...]:
    """Return the innermost active procedure for every physical line."""
    procedures = sorted(procedure_cases, key=lambda procedure: procedure.start_line)
    active: list[ProcedureDeclarationCases] = []
    result: list[ProcedureDeclarationCases | None] = []
    next_procedure = 0
    for line_number in range(line_count):
        while active and active[-1].end_line < line_number:
            active.pop()
        while next_procedure < len(procedures) and procedures[next_procedure].start_line <= line_number:
            active.append(procedures[next_procedure])
            next_procedure += 1
        while active and active[-1].end_line < line_number:
            active.pop()
        result.append(active[-1] if active else None)
    return tuple(result)


def _declared_names_by_line(scoped_names: Iterable[ScopedDeclaredNames], line_count: int) -> tuple[frozenset[str], ...]:
    """Return the union of declared names visible on every physical line."""
    starts: dict[int, list[frozenset[str]]] = {}
    ends: dict[int, list[frozenset[str]]] = {}
    for scope in scoped_names:
        starts.setdefault(scope.start_line, []).append(scope.names)
        ends.setdefault(scope.end_line, []).append(scope.names)

    counts: dict[str, int] = {}
    result: list[frozenset[str]] = []
    for line_number in range(line_count):
        for names in starts.get(line_number, ()):
            for name in names:
                counts[name] = counts.get(name, 0) + 1
        result.append(frozenset(counts))
        for names in ends.get(line_number, ()):
            for name in names:
                remaining = counts[name] - 1
                if remaining:
                    counts[name] = remaining
                else:
                    del counts[name]
    return tuple(result)


def extract_named_scope_cases(source: str) -> tuple[NamedScopeCase, ...]:
    """Extract start-name cases for named module, program, and procedure scopes."""
    records: list[dict[str, object]] = []
    stack: list[int] = []
    statements = _code_statements(source)
    last_line = statements[-1].end_line if statements else 0
    for statement in statements:
        text = statement.text
        end = NAMED_SCOPE_END.match(text)
        if end or BARE_PROGRAM_UNIT_END.match(text):
            kind = end.group(1).lower() if end else None
            for index in range(len(stack) - 1, -1, -1):
                record = records[stack[index]]
                is_bare_program_unit_end = kind is None and record["kind"] in {"function", "subroutine", "program"}
                if is_bare_program_unit_end or record["kind"] == kind:
                    record["end_line"] = statement.end_line
                    del stack[index]
                    break
            continue

        header = _scope_header(text)
        module = MODULE_DECLARATION.match(text)
        if header:
            kind, name = header[0], header[1]
        elif module:
            kind, name = "module", module.group(1)
        else:
            continue
        records.append({"kind": kind, "name": name, "start_line": statement.start_line, "end_line": last_line})
        stack.append(len(records) - 1)

    return tuple(
        NamedScopeCase(record["kind"], int(record["start_line"]), int(record["end_line"]), record["name"])
        for record in records
    )


def _resolve_case_occurrences(occurrences: Mapping[str, list[str]]) -> dict[str, str]:
    """Return names whose declarations all use one exact spelling."""
    return {normalized: spellings[0] for normalized, spellings in occurrences.items() if len(set(spellings)) == 1}


def _resolve_type_occurrences(occurrences: Mapping[str, set[str]]) -> dict[str, str]:
    """Return type names whose declarations agree case-insensitively."""
    return {
        normalized: next(iter(type_names))
        for normalized, type_names in occurrences.items()
        if len({type_name.lower() for type_name in type_names}) == 1
    }


def extract_procedure_cases(source: str) -> tuple[ProcedureDeclarationCases, ...]:
    """Extract local variable cases for each function or subroutine in *source*."""
    records: list[dict[str, object]] = []
    stack: list[int] = []
    statements = _code_statements(source)
    for statement in statements:
        text = statement.text
        if PROCEDURE_UNIT_END.match(text):
            if stack:
                records[stack.pop()]["end_line"] = statement.end_line
            continue
        if BARE_PROGRAM_UNIT_END.match(text):
            if stack:
                records[stack.pop()]["end_line"] = statement.end_line
            continue

        header = _scope_header(text)
        if header:
            _, procedure_name, header_names = header
            records.append(
                {
                    "start_line": statement.start_line,
                    "end_line": statements[-1].end_line if statements else statement.end_line,
                    "header_names": header_names,
                    "declarations": [],
                    "types": {},
                    "contains": False,
                    "name": procedure_name,
                }
            )
            stack.append(len(records) - 1)
            continue

        if not stack:
            continue
        record = records[stack[-1]]
        if CONTAINS_STATEMENT.match(text):
            record["contains"] = True
        elif not record["contains"]:
            record["declarations"].extend(_declared_variable_names(text))
            record["types"].update(_declared_variable_types(text))
            alias = _select_type_alias(text)
            if alias is not None:
                record["declarations"].append(alias)

    local_cases: list[ProcedureDeclarationCases] = []
    for record in records:
        explicit_names = record["declarations"]
        header_names = record["header_names"]
        occurrences: dict[str, list[str]] = {}
        for name in explicit_names:
            occurrences.setdefault(name.lower(), []).append(name)
        explicit_keys = set(occurrences)
        for name in header_names:
            if name.lower() not in explicit_keys:
                occurrences.setdefault(name.lower(), []).append(name)
        local_cases.append(
            ProcedureDeclarationCases(
                int(record["start_line"]),
                int(record["end_line"]),
                _resolve_case_occurrences(occurrences),
                frozenset(occurrences),
                record["types"],
            )
        )
    return tuple(local_cases)


@dataclass(frozen=True)
class DeclarationSummary:
    """Project-wide declaration facts extracted in one logical-statement pass."""

    declared_names: tuple[DeclaredName, ...]
    module_variable_names: tuple[str, ...]
    module_variable_types: tuple[tuple[str, str], ...]
    type_bound_procedure_names: tuple[str, ...]
    type_component_names: tuple[tuple[str, str], ...]
    type_component_types: tuple[tuple[str, str, str], ...]


@lru_cache(maxsize=128)
def _declaration_summary(source: str) -> DeclarationSummary:
    """Extract declaration facts once instead of rescanning the same statements."""
    declarations: list[DeclaredName] = []
    module_variable_names: list[str] = []
    module_variable_types: dict[str, str] = {}
    type_bound_procedure_names: list[str] = []
    type_component_names: list[tuple[str, str]] = []
    type_component_types: list[tuple[str, str, str]] = []

    in_module = False
    module_contains = False
    module_type_depth = 0
    module_interface_depth = 0
    type_stack: list[str] = []

    for statement in _code_statements(source):
        text = statement.text
        module = MODULE_DECLARATION.match(text)
        type_name = _declared_type_name(text)
        type_start = type_name is not None
        type_end = TYPE_DEFINITION_END.match(text) is not None
        procedure_declaration = DECLARATION_PROCEDURE_START.match(text) is not None
        variable_names = _declared_variable_names(text)

        # General module/type/procedure declarations.
        if END_STATEMENT.match(text) is None:
            if module:
                declarations.append(DeclaredName("module", module.group(1)))
            if type_name:
                declarations.append(DeclaredName("type", type_name))
            procedure = PROCEDURE_DECLARATION_NAME.search(text)
            if procedure:
                declarations.append(DeclaredName("procedure", procedure.group(2)))

        # Module specification-part names and types. This intentionally mirrors
        # the previous extractor's treatment of ordinary derived-type components
        # as globally available symbol-case evidence while excluding bindings.
        if module:
            in_module = True
            module_contains = False
            module_type_depth = 0
            module_interface_depth = 0
        elif in_module and END_MODULE.match(text):
            in_module = False
        elif in_module and not module_contains:
            if INTERFACE_END.match(text):
                module_interface_depth = max(0, module_interface_depth - 1)
            elif INTERFACE_START.match(text):
                module_interface_depth += 1
            elif module_interface_depth:
                pass
            elif type_start:
                module_type_depth += 1
            elif module_type_depth:
                if type_end:
                    module_type_depth -= 1
                elif not procedure_declaration:
                    module_variable_names.extend(variable_names)
            elif CONTAINS_STATEMENT.match(text):
                module_contains = True
            else:
                module_variable_names.extend(variable_names)
                module_variable_types.update(_declared_variable_types(text))

        # Type-bound procedure and component facts for every derived-type body.
        if type_start:
            type_stack.append(type_name)
        elif type_stack:
            if type_end:
                type_stack.pop()
            elif procedure_declaration:
                type_bound_procedure_names.extend(variable_names)
            else:
                owner_type = type_stack[-1]
                type_component_names.extend((owner_type, name) for name in variable_names)
                type_component_types.extend(
                    (owner_type, name, component_type)
                    for name, component_type in _declared_variable_types(text).items()
                )

    return DeclarationSummary(
        tuple(declarations),
        tuple(module_variable_names),
        tuple(module_variable_types.items()),
        tuple(type_bound_procedure_names),
        tuple(type_component_names),
        tuple(type_component_types),
    )


def extract_module_variable_names(source: str) -> list[str]:
    """Extract variable names declared in module specification parts."""
    return list(_declaration_summary(source).module_variable_names)


def extract_module_variable_types(source: str) -> dict[str, str]:
    """Extract derived types for variables in module specification parts."""
    return dict(_declaration_summary(source).module_variable_types)


def extract_type_bound_procedure_names(source: str) -> list[str]:
    """Extract type-bound procedure binding names from module type definitions."""
    return list(_declaration_summary(source).type_bound_procedure_names)


def extract_type_component_names(source: str) -> list[tuple[str, str]]:
    """Extract ordinary derived-type component names with their type."""
    return list(_declaration_summary(source).type_component_names)


def extract_type_component_types(source: str) -> list[tuple[str, str, str]]:
    """Extract derived-type component types for resolving chained members."""
    return list(_declaration_summary(source).type_component_types)


def extract_declared_names(source: str) -> list[DeclaredName]:
    """Extract module, derived-type, function, and subroutine declaration names."""
    return list(_declaration_summary(source).declared_names)


def extract_scoped_declared_names(source: str) -> tuple[ScopedDeclaredNames, ...]:
    """Extract type, procedure, and module-variable names for each enclosing
    module, program, function, or subroutine in *source*.

    A name is attributed to the scope that textually encloses its
    declaration -- not to the scope it may itself open, for a procedure
    declared by its own header line -- so that a private module variable, a
    module-contained procedure, or a locally-defined type is visible only
    where Fortran actually allows it to be referenced. Declarations that sit
    outside every module/program/function/subroutine (rare: a legacy
    top-level file) are dropped rather than treated as file-wide, since
    ``local_names`` already covers ordinary procedure-local declarations and
    there is no broader scope left to attribute them to safely.
    """
    records: list[dict[str, object]] = []
    stack: list[int] = []
    statements = _code_statements(source)
    last_line = statements[-1].end_line if statements else 0

    def add_name(name: str) -> None:
        if stack:
            records[stack[-1]]["names"].add(name.lower())

    for statement in statements:
        text = statement.text
        end = NAMED_SCOPE_END.match(text)
        if end or BARE_PROGRAM_UNIT_END.match(text):
            kind = end.group(1).lower() if end else None
            for index in range(len(stack) - 1, -1, -1):
                record = records[stack[index]]
                is_bare_program_unit_end = kind is None and record["kind"] in {"function", "subroutine", "program"}
                if is_bare_program_unit_end or record["kind"] == kind:
                    record["end_line"] = statement.end_line
                    del stack[index]
                    break
            continue

        header = _scope_header(text)
        module = MODULE_DECLARATION.match(text)
        if header:
            kind, name = header[0], header[1]
        elif module:
            kind, name = "module", module.group(1)
        else:
            kind = name = None

        if kind is not None:
            # Attributed to the *enclosing* scope (or dropped, if there is
            # none): a procedure is not a member of the scope it opens.
            add_name(name)
            records.append(
                {
                    "kind": kind,
                    "start_line": statement.start_line,
                    "end_line": last_line,
                    "names": set(),
                    "contains_seen": False,
                    "type_depth": 0,
                    "interface_depth": 0,
                }
            )
            stack.append(len(records) - 1)
            continue

        if not stack:
            continue
        frame = records[stack[-1]]
        if INTERFACE_END.match(text):
            frame["interface_depth"] = max(0, frame["interface_depth"] - 1)
            continue
        if INTERFACE_START.match(text):
            frame["interface_depth"] += 1
            continue
        if frame["interface_depth"]:
            continue
        type_name = _declared_type_name(text)
        if TYPE_DECLARATION_START.match(text) and type_name:
            add_name(type_name)
            frame["type_depth"] += 1
            continue
        if frame["type_depth"]:
            if TYPE_DEFINITION_END.match(text):
                frame["type_depth"] -= 1
            # A type body's own `contains` (introducing type-bound procedures)
            # must not be mistaken for the enclosing scope's `contains`, so
            # this check has to stay ahead of the `contains` check below.
            # Components require `%` access and are handled separately; they
            # are not candidates for overriding a bare token's lowercasing.
            continue
        if CONTAINS_STATEMENT.match(text):
            frame["contains_seen"] = True
            continue
        if frame["kind"] == "module" and not frame["contains_seen"]:
            for name in _declared_variable_names(text):
                add_name(name)

    return tuple(
        ScopedDeclaredNames(int(record["start_line"]), int(record["end_line"]), frozenset(record["names"]))
        for record in records
    )


def declared_names_at(scoped_names: Iterable[ScopedDeclaredNames], line_number: int) -> frozenset[str]:
    """Return names in scope at *line_number*, from every enclosing scope.

    A name declared at module level stays visible in every procedure nested
    inside that module, so this unions every scope whose range contains
    *line_number* rather than returning only the innermost one.
    """
    names: set[str] = set()
    for scope in scoped_names:
        if scope.start_line <= line_number <= scope.end_line:
            names |= scope.names
    return frozenset(names)


def _case_for_file(
    path: Path,
    declarations: Mapping[tuple[str, str], list[tuple[Path, str]]],
    kind: str,
) -> dict[str, str]:
    """Choose unambiguous declaration spellings for *path* and *kind*."""
    selected: dict[str, str] = {}
    for (declaration_kind, normalized), occurrences in declarations.items():
        if declaration_kind != kind:
            continue
        local = [spelling for declaration_path, spelling in occurrences if declaration_path == path]
        local_cases = set(local)
        all_cases = {spelling for _, spelling in occurrences}
        if local_cases and len(local_cases) == 1:
            selected[normalized] = local[0]
        elif not local and len(all_cases) == 1:
            selected[normalized] = occurrences[0][1]
    return selected


def collect_declaration_cases(
    sources: Mapping[Path, str],
    target_paths: Collection[Path] | None = None,
) -> dict[Path, FileDeclarationCases]:
    """Collect project declarations and resolve casing for requested target files.

    All *sources* participate in project-wide declaration resolution, but local
    procedure/scope metadata is only needed for files that will actually be
    formatted. Omitting *target_paths* preserves the original all-files API.

    A declaration in the current file wins over declarations in other files. A
    name declared with different spellings more than once in the current file,
    or ambiguously elsewhere with no current-file declaration, is omitted.
    """
    targets = tuple(sources) if target_paths is None else tuple(dict.fromkeys(target_paths))
    missing = [path for path in targets if path not in sources]
    if missing:
        raise KeyError(f"Target paths are missing from sources: {missing}")

    declarations: dict[tuple[str, str], list[tuple[Path, str]]] = {}
    variable_type_occurrences: dict[str, set[str]] = {}
    component_type_occurrences: dict[tuple[str, str], set[str]] = {}
    procedure_cases = {path: extract_procedure_cases(sources[path]) for path in targets}
    scope_cases = {path: extract_named_scope_cases(sources[path]) for path in targets}
    for path, source in sources.items():
        summary = _declaration_summary(source)
        for declaration in summary.declared_names:
            kind = "module" if declaration.kind == "module" else "symbol"
            key = (kind, declaration.name.lower())
            declarations.setdefault(key, []).append((path, declaration.name))
        for name in summary.module_variable_names:
            key = ("symbol", name.lower())
            declarations.setdefault(key, []).append((path, name))
        for name in summary.type_bound_procedure_names:
            key = ("type_procedure", name.lower())
            declarations.setdefault(key, []).append((path, name))
        for type_name, name in summary.type_component_names:
            key = ("type_component", f"{type_name.lower()}\0{name.lower()}")
            declarations.setdefault(key, []).append((path, name))
        for name, type_name in summary.module_variable_types:
            variable_type_occurrences.setdefault(name, set()).add(type_name)
        for type_name, name, component_type in summary.type_component_types:
            component_type_occurrences.setdefault((type_name.lower(), name.lower()), set()).add(component_type)

    variable_type_cases = _resolve_type_occurrences(variable_type_occurrences)
    type_component_type_cases = _resolve_type_occurrences(
        {f"{type_name}\0{name}": type_names for (type_name, name), type_names in component_type_occurrences.items()}
    )

    cases: dict[Path, FileDeclarationCases] = {}
    for path in targets:
        module_cases = _case_for_file(path, declarations, "module")
        # Procedure-local names are handled while replacing tokens, where their
        # scope is known. Filtering them here would also suppress a global
        # component spelling everywhere else in this file (e.g. P%H0 when a
        # different procedure has a local H0 argument).
        symbol_cases = _case_for_file(path, declarations, "symbol")
        type_procedure_cases = _case_for_file(path, declarations, "type_procedure")
        type_component_cases = {
            tuple(key.split("\0", 1)): spelling
            for key, spelling in _case_for_file(path, declarations, "type_component").items()
        }
        resolved_component_types = {
            tuple(key.split("\0", 1)): type_name for key, type_name in type_component_type_cases.items()
        }
        cases[path] = FileDeclarationCases(
            module_cases,
            symbol_cases,
            procedure_cases[path],
            scope_cases[path],
            type_procedure_cases,
            type_component_cases,
            variable_type_cases,
            resolved_component_types,
        )
    return cases


def member_owner_type(
    context: str,
    index: int,
    local_types: Mapping[str, str],
    variable_types: Mapping[str, str],
    type_component_types: Mapping[tuple[str, str], str],
) -> str | None:
    """Resolve the type of the member expression immediately before a component."""
    position = index - 1
    while position >= 0 and context[position].isspace():
        position -= 1
    if position < 0 or context[position] != "%":
        return None
    position -= 1

    reversed_chain: list[str] = []
    while position >= 0:
        while position >= 0 and context[position].isspace():
            position -= 1
        end = position + 1
        while position >= 0 and _is_identifier_char(context[position]):
            position -= 1
        start = position + 1
        if start == end or not _is_identifier_start(context[start]):
            return None
        reversed_chain.append(context[start:end])
        while position >= 0 and context[position].isspace():
            position -= 1
        if position < 0 or context[position] != "%":
            break
        position -= 1

    chain = reversed(reversed_chain)
    first = next(chain, None)
    if first is None:
        return None
    owner_type = local_types.get(first.lower()) or variable_types.get(first.lower())
    for component in chain:
        if owner_type is None:
            return None
        owner_type = type_component_types.get((owner_type.lower(), component.lower()))
    return owner_type


def replace_declared_cases(
    source: str,
    module_cases: Mapping[str, str],
    symbol_cases: Mapping[str, str],
    procedure_cases: Iterable[ProcedureDeclarationCases] = (),
    scope_cases: Iterable[NamedScopeCase] = (),
    type_procedure_cases: Mapping[str, str] | None = None,
    type_component_cases: Mapping[tuple[str, str], str] | None = None,
    variable_type_cases: Mapping[str, str] | None = None,
    type_component_type_cases: Mapping[tuple[str, str], str] | None = None,
    preprocessor_names: Collection[str] = frozenset(),
) -> str:
    """Match declaration case in USE statements and references to declared symbols."""
    procedure_cases = tuple(procedure_cases)
    scope_cases = tuple(scope_cases)
    variable_types = extract_variable_types(source)
    variable_type_cases = variable_type_cases or {}
    resolved_variable_types = {**variable_type_cases, **variable_types}
    type_procedure_cases = type_procedure_cases or {}
    type_component_cases = type_component_cases or {}
    type_component_type_cases = type_component_type_cases or {}
    if (
        not module_cases
        and not symbol_cases
        and not procedure_cases
        and not scope_cases
        and not type_procedure_cases
        and not type_component_cases
        and not variable_type_cases
        and not type_component_type_cases
    ):
        return source

    output: list[str] = []
    prefix = ""
    _, states = _scan_source(source)
    active_procedures = _active_procedures_by_line(procedure_cases, len(states))
    end_scope_by_line = {scope.end_line: scope for scope in scope_cases}
    for state in states:
        line_number = state.number
        line = state.text
        if state.is_cpp:
            output.append(line)
            prefix = ""
            continue

        active_procedure = active_procedures[line_number]
        local_cases = active_procedure.local_cases if active_procedure else {}
        local_names = active_procedure.local_names if active_procedure else frozenset()
        starting_quote = state.quote_in
        quote = starting_quote
        context_line = blank_leading_continuation(line) if prefix else line
        context = prefix + context_line
        masked_context = code_context(context)
        use_match = USE_MODULE.match(context)
        end_match = NAMED_SCOPE_END.match(context)
        module_start = use_match.start(1) if use_match else -1
        module_end = use_match.end(1) if use_match else -1
        end_scope = end_scope_by_line.get(line_number)
        end_name_start = end_match.start(2) if end_match and end_match.group(2) else -1
        end_name_end = end_match.end(2) if end_match and end_match.group(2) else -1
        line_output: list[str] = []
        index = 0
        while index < len(line):
            char = line[index]
            if quote:
                line_output.append(char)
                if char == quote:
                    if line[index + 1 : index + 2] == quote:
                        line_output.append(line[index + 1])
                        index += 2
                        continue
                    quote = None
                index += 1
            elif char in "\"'":
                quote = char
                line_output.append(char)
                index += 1
            elif char == "!":
                line_output.append(line[index:])
                break
            elif _is_identifier_start(char):
                token_end = _identifier_end(line, index)
                token = line[index:token_end]
                normalized = token.lower()
                absolute_start = len(prefix) + index
                member_component = _preceded_by_percent(context, absolute_start)
                if (
                    (
                        not member_component
                        and masked_context.rfind("::", 0, absolute_start) >= 0
                        and is_contextual_identifier(context, absolute_start, context=masked_context)
                    )
                    or token in preprocessor_names
                    or (normalized in {"e", "d"} and is_real_literal_exponent_marker(context, absolute_start, token))
                ):
                    replacement = None
                elif (
                    end_scope
                    and end_match
                    and end_scope.kind == end_match.group(1).lower()
                    and end_name_start <= absolute_start < end_name_end
                ):
                    replacement = end_scope.name
                elif member_component:
                    owner_type = member_owner_type(
                        context,
                        absolute_start,
                        local_types=active_procedure.local_types if active_procedure else {},
                        variable_types=resolved_variable_types,
                        type_component_types=type_component_type_cases,
                    )
                    replacement = (
                        type_procedure_cases.get(normalized)
                        or (type_component_cases.get((owner_type.lower(), normalized)) if owner_type else None)
                        or symbol_cases.get(normalized)
                    )
                elif module_start <= absolute_start < module_end:
                    replacement = module_cases.get(normalized)
                elif use_match:
                    replacement = symbol_cases.get(normalized)
                else:
                    replacement = local_cases.get(normalized)
                    if replacement is None and normalized not in local_names and normalized not in INTRINSIC_NAMES:
                        replacement = symbol_cases.get(normalized)
                line_output.append(replacement or token)
                index = token_end
            else:
                line_output.append(char)
                index += 1
        updated_line = "".join(line_output)
        output.append(updated_line)
        if not state.is_blank and not state.is_comment:
            prefix = statement_context(prefix, line, starting_quote) if state.continuation_out else ""
    return "".join(output)


def _is_identifier_start(char: str) -> bool:
    return "A" <= char <= "Z" or "a" <= char <= "z"


def _is_identifier_char(char: str) -> bool:
    return _is_identifier_start(char) or "0" <= char <= "9" or char == "_"


def _identifier_end(text: str, index: int) -> int:
    end = index + 1
    length = len(text)
    while end < length:
        char = text[end]
        if not ("A" <= char <= "Z" or "a" <= char <= "z" or "0" <= char <= "9" or char == "_"):
            break
        end += 1
    return end


def _preceded_by_percent(text: str, index: int) -> bool:
    index -= 1
    while index >= 0 and text[index].isspace():
        index -= 1
    return index >= 0 and text[index] == "%"


def is_real_literal_exponent_marker(context: str, index: int, token: str) -> bool:
    """Return whether the token is the exponent marker in a real literal."""
    if token.lower() not in {"e", "d"}:
        return False
    mantissa = context[:index]
    exponent = context[index + len(token) :]
    return re.search(r"(?:\d+(?:\.\d*)?|\.\d+)$", mantissa) is not None and re.match(r"[+-]?\d", exponent) is not None


@lru_cache(maxsize=8192)
def code_context(line: str) -> str:
    """Return code from *line* with quoted text and comments masked."""
    if "'" not in line and '"' not in line:
        comment_start = line.find("!")
        return line if comment_start < 0 else line[:comment_start]
    context = mask_quoted_text(line)
    comment_start = context.find("!")
    return context if comment_start < 0 else context[:comment_start]


def is_contextual_identifier(line: str, index: int, *, context: str | None = None) -> bool:
    """Return whether the token at *index* is being used as an identifier."""
    context = code_context(line) if context is None else context

    previous = index - 1
    while previous >= 0 and context[previous].isspace():
        previous -= 1
    if previous >= 0 and context[previous] == "%":
        return True

    declaration_start = context.rfind("::", 0, index)
    if declaration_start < 0:
        return False
    item_start = declaration_start + 2
    depth = 0
    for position in range(item_start, index):
        char = context[position]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            item_start = position + 1

    # An initializer/association makes the token a value expression rather than
    # the name being declared.
    for position in range(item_start, index):
        if context[position] != "=":
            continue
        previous_char = context[position - 1] if position > item_start else ""
        following = context[position + 1 : position + 2]
        if following == ">" or (previous_char not in "<>=/" and following not in {"=", ">"}):
            return False
    return True


def lowercase_keyword(
    token: str,
    line: str,
    index: int,
    preserved_names: Collection[str] = frozenset(),
    local_names: Collection[str] = frozenset(),
    preprocessor_names: Collection[str] = frozenset(),
    file_declared_names: Collection[str] = frozenset(),
    uppercase_single_l: bool = False,
    *,
    context: str | None = None,
) -> str:
    """Lowercase language and intrinsic names without overriding declarations."""
    token_lower = token.lower()
    context = code_context(line) if context is None else context
    token_end = index + len(token)

    if token_lower in {"e", "d"} and is_real_literal_exponent_marker(context, index, token):
        return token_lower
    if _preceded_by_percent(context, index):
        return token
    if context.rfind("::", 0, index) >= 0 and is_contextual_identifier(line, index, context=context):
        return token
    if token in preprocessor_names:
        return token

    if token_lower in local_names:
        return token
    if token_lower in file_declared_names and not is_specifier_keyword_argument(context, token_end):
        return token

    def next_nonspace(start: int) -> int:
        while start < len(context) and context[start] in " \t":
            start += 1
        return start

    def previous_word_is_double() -> bool:
        end = index
        while end > 0 and context[end - 1].isspace():
            end -= 1
        start = end
        while start > 0 and _is_identifier_char(context[start - 1]):
            start -= 1
        return context[start:end].lower() == "double"

    if token_lower in FORTRAN_STANDARD_WORDS:
        if token_lower in DECLARATION_ATTRIBUTES and context.find("::", token_end) < 0:
            return token
        if token_lower == "only":
            following = next_nonspace(token_end)
            if following >= len(context) or context[following] != ":":
                return token
        elif token_lower == "bind":
            following = next_nonspace(token_end)
            if following >= len(context) or context[following] != "(":
                return token
            following = next_nonspace(following + 1)
            if following >= len(context) or context[following].lower() != "c":
                return token
            following = next_nonspace(following + 1)
            if following >= len(context) or context[following] != ")":
                return token
        elif token_lower == "kind":
            following = next_nonspace(token_end)
            if following >= len(context) or context[following] not in "(=":
                return token
        elif token_lower == "precision" and not previous_word_is_double():
            return token
        return token_lower

    if token_lower in PARENTHESIZED_STATEMENT_NAMES and token_lower not in local_names:
        following = next_nonspace(token_end)
        if following < len(context) and context[following] == "(":
            return token_lower

    if uppercase_single_l and token_lower == "l":
        return "L"
    if token_lower in INTRINSIC_NAMES:
        if token_lower == "precision":
            following = next_nonspace(token_end)
            if (following >= len(context) or context[following] != "(") and not previous_word_is_double():
                return token
        return token_lower
    if token_lower in preserved_names:
        return token
    return token


def is_binary_arithmetic_operator(line: str, index: int, operator: str) -> bool:
    """Return whether an arithmetic token is infix rather than a unary sign."""
    previous = index - 1
    while previous >= 0 and line[previous].isspace():
        previous -= 1
    following = index + len(operator)
    while following < len(line) and line[following].isspace():
        following += 1

    if operator == "//":
        return following < len(line)
    if previous < 0 or following >= len(line):
        return False
    if line[previous] not in ")]_.0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
        return False
    # A trailing `.` closes a dotted operator (`x.lt.-1`) as often as it ends a
    # real literal (`1.-x`); only the latter leaves an operand to the left.
    if line[previous] == "." and re.search(r"\.[A-Za-z]+\.$", line[: previous + 1]):
        return False
    if operator in "+-" and re.search(r"(?:\d|\.)[eEdD]$", line[:index]):
        return False
    return line[following] not in "),]"


def is_named_parameter(line: str, index: int) -> bool:
    """Return whether ``=`` at *index* is a named argument or specification."""
    end = index
    while end > 0 and line[end - 1].isspace():
        end -= 1
    start = end
    while start > 0 and _is_identifier_char(line[start - 1]):
        start -= 1
    if start == end or not _is_identifier_start(line[start]):
        return False
    prefix = start
    while prefix > 0 and line[prefix - 1].isspace():
        prefix -= 1
    if prefix == 0 or line[prefix - 1] not in "(,":
        return False

    depth = 0
    quote: str | None = None
    position = 0
    while position < index:
        char = line[position]
        if quote:
            if char == quote:
                if line[position + 1 : position + 2] == quote:
                    position += 2
                    continue
                quote = None
        elif char in "\"'":
            quote = char
        elif char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
        position += 1
    return depth > 0


def is_specifier_keyword_argument(context: str, token_end: int) -> bool:
    """Return whether the token ending at *token_end* names a fixed ``keyword=`` specifier."""
    equals_index = token_end
    while equals_index < len(context) and context[equals_index] in " \t":
        equals_index += 1
    if equals_index >= len(context) or context[equals_index] != "=":
        return False
    following = context[equals_index + 1 : equals_index + 2]
    if following in {"=", ">"}:
        return False
    return is_named_parameter(context, equals_index)


def pop_trailing_whitespace(output: list[str]) -> str:
    """Remove and return the whitespace at the end of *output*."""
    trailing: list[str] = []
    while output and output[-1].isspace():
        trailing.append(output.pop())
    trailing.reverse()
    return "".join(trailing)


def append_spaced_operator(output: list[str], operator: str) -> None:
    """Append *operator* with one space on either side, keeping any line indent.

    A continuation line may start with an operator (``.and.``, ``==``, ...), in
    which case the whitespace before it is the line's indentation rather than
    operator padding and must survive unchanged.
    """
    indent = pop_trailing_whitespace(output)
    output.append(" " if output else indent)
    output.extend((operator, " "))


def append_compact_operator(output: list[str], operator: str) -> None:
    """Append *operator* with no surrounding whitespace, keeping any line indent."""
    indent = pop_trailing_whitespace(output)
    if not output:
        output.append(indent)
    output.append(operator)


def skip_operator_whitespace(line: str, index: int) -> int:
    """Skip spaces after an operator without consuming a line ending."""
    while index < len(line) and line[index].isspace() and line[index] != "\n":
        index += 1
    return index


def append_normalized_operator(output: list[str], line: str, index: int, context: str, offset: int) -> int | None:
    """Append the normalized operator at *index*, if there is one."""
    char = line[index]
    if char not in ".=<>+-*/":
        return None

    if char == ".":
        operator = LEGACY_OPERATOR.match(line, index)
        if operator:
            append_spaced_operator(output, MODERN_OPERATOR[operator.group(1).lower()])
            return skip_operator_whitespace(line, operator.end())
        operator = SPACED_OPERATOR.match(line, index)
        if operator:
            operator_text = operator.group().lower()
            append_spaced_operator(output, operator_text)
            return skip_operator_whitespace(line, operator.end())
        return None

    # Relational/assignment operators overlap arithmetic '/' so try this tier
    # first only for characters that can actually begin one.
    if char in "=<>/":
        operator = SPACED_OPERATOR.match(line, index)
        if operator:
            operator_text = operator.group()
            if operator_text == "=" and is_named_parameter(context, offset + index):
                append_compact_operator(output, operator_text)
            else:
                append_spaced_operator(output, operator_text)
            return skip_operator_whitespace(line, operator.end())

    operator = ARITHMETIC_OPERATOR.match(line, index)
    if operator:
        operator_text = operator.group()
        if not is_binary_arithmetic_operator(context, offset + index, operator_text):
            output.append(operator_text)
            return skip_operator_whitespace(line, operator.end())
        if operator_text in COMPACT_ARITHMETIC_OPERATORS:
            append_compact_operator(output, operator_text)
        else:
            append_spaced_operator(output, operator_text)
        return skip_operator_whitespace(line, operator.end())
    return None


def uppercase_openmp_keywords(line: str, preprocessor_names: Collection[str] = frozenset()) -> str:
    """Uppercase OpenMP vocabulary without changing case-sensitive macro names."""
    macro_names = {name.lower() for name in preprocessor_names}
    output: list[str] = []
    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote:
            output.append(char)
            if char == quote:
                if line[index + 1 : index + 2] == quote:
                    output.append(line[index + 1])
                    index += 2
                    continue
                quote = None
            index += 1
        elif char in "\"'":
            quote = char
            output.append(char)
            index += 1
        elif char == "!":
            output.append(line[index:])
            break
        elif _is_identifier_start(char):
            end = _identifier_end(line, index)
            token = line[index:end]
            normalized = token.lower()
            replacement = (
                token if normalized in macro_names else token.upper() if normalized in OPENMP_KEYWORDS else token
            )
            output.append(replacement)
            index = end
        else:
            output.append(char)
            index += 1
    return "".join(output)


def normalize_openmp_clause_separators(line: str) -> str:
    """Separate adjacent OpenMP clauses without touching literals or comments."""
    output: list[str] = []
    quote: str | None = None
    index = 0
    length = len(line)
    while index < length:
        char = line[index]
        if quote:
            output.append(char)
            if char == quote:
                if line[index + 1 : index + 2] == quote:
                    output.append(line[index + 1])
                    index += 2
                    continue
                quote = None
            index += 1
            continue
        if char in "\"'":
            quote = char
            output.append(char)
            index += 1
            continue
        if char == "!":
            output.append(line[index:])
            break
        if char == ")":
            whitespace_end = index + 1
            while whitespace_end < length and line[whitespace_end] in " \t":
                whitespace_end += 1
            token_end = whitespace_end
            if token_end < length and _is_identifier_start(line[token_end]):
                token_end = _identifier_end(line, token_end)
                opening = token_end
                while opening < length and line[opening] in " \t":
                    opening += 1
                if opening < length and line[opening] == "(" and whitespace_end > index + 1:
                    output.append("), ")
                    index = whitespace_end
                    continue
        output.append(char)
        index += 1
    return "".join(output)


def lowercase_line(
    line: str,
    quote: str | None = None,
    prefix: str = "",
    preserved_names: Collection[str] = frozenset(),
    local_names: Collection[str] = frozenset(),
    preprocessor_names: Collection[str] = frozenset(),
    file_declared_names: Collection[str] = frozenset(),
    uppercase_single_l: bool = False,
) -> tuple[str, str | None]:
    """Lowercase keywords and normalize code operators outside quoted text."""
    if is_preprocessor_line(line):
        return line, quote

    stripped = line.lstrip(" \t")
    if stripped.startswith("!$"):
        sentinel_end = len(line) - len(stripped) + 2
        sentinel = line[:sentinel_end]
        directive = line[sentinel_end:]
        directive_body = directive.lstrip(" \t")
        if directive_body[:3].lower() == "omp" and (
            len(directive_body) == 3 or not _is_identifier_char(directive_body[3])
        ):
            normalized = uppercase_openmp_keywords(directive, preprocessor_names)
            normalized, _ = normalize_delimiter_spacing(normalized)
            return sentinel + normalize_openmp_clause_separators(normalized), quote
        normalized, quote = lowercase_line(
            directive,
            quote,
            preserved_names=preserved_names,
            local_names=local_names,
            preprocessor_names=preprocessor_names,
            file_declared_names=file_declared_names,
            uppercase_single_l=uppercase_single_l,
        )
        normalized, quote = normalize_keyword_spacing(normalized, quote)
        normalized, _ = normalize_delimiter_spacing(normalized)
        return sentinel + normalized, quote

    context = prefix + (blank_leading_continuation(line) if prefix else line)
    masked_context = code_context(context)
    offset = len(prefix)
    output: list[str] = []
    index = 0
    while index < len(line):
        char = line[index]
        if quote:
            output.append(char)
            if char == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    output.append(line[index + 1])
                    index += 2
                    continue
                quote = None
            index += 1
            continue

        if char in "\"'":
            quote = char
            output.append(char)
            index += 1
            continue
        if char == "!":
            comment = line[index:]
            output.append(format_comment_operators(comment) if comment_contains_assignment(comment) else comment)
            break

        # Numeric exponent matching is meaningful only at a numeric literal
        # start; the previous implementation attempted this regex at almost
        # every character in every line.
        if "0" <= char <= "9" or char == ".":
            real_literal = REAL_LITERAL_EXPONENT.match(line, index)
            if real_literal:
                output.append(
                    real_literal.group("mantissa")
                    + real_literal.group("marker").lower()
                    + real_literal.group("exponent")
                )
                index = real_literal.end()
                continue

        if char in ".=<>+-*/":
            operator_end = append_normalized_operator(output, line, index, context, offset)
            if operator_end is not None:
                index = operator_end
                continue

        if _is_identifier_start(char):
            end = _identifier_end(line, index)
            token = line[index:end]
            output.append(
                lowercase_keyword(
                    token,
                    context,
                    offset + index,
                    preserved_names,
                    local_names,
                    preprocessor_names,
                    file_declared_names=file_declared_names,
                    uppercase_single_l=uppercase_single_l,
                    context=masked_context,
                )
            )
            index = end
            continue

        output.append(char)
        index += 1
    return "".join(output), quote


def _modernize_array_constructor_delimiters(code: str) -> str:
    """Modernize old ``(/ ... /)`` array constructors, but never FORMAT syntax."""
    # Slash is an edit descriptor inside FORMAT(...), where the same character
    # sequence is legal and has completely different semantics.
    if FORMAT_STATEMENT_START.match(code):
        return code
    code = re.sub(r"\(\s*/\s*", "[", code)
    return re.sub(r"/\s*\)", "]", code)


def _normalize_keyword_spacing_code(code: str) -> str:
    """Normalize keyword/layout spacing in an unquoted code segment."""
    if not code:
        return code

    lowered = code.lower()

    # Most lines need only a tiny subset of these rules. Cheap substring guards
    # avoid dispatching dozens of regex substitutions for unrelated statements.
    if "common" in lowered and "/" in code:
        code = COMMON_BLOCK_PREFIX.sub(lambda match: f"{match.group(1)} /{match.group(2)}/ ", code, count=1)
    if "(" in code and "/" in code:
        code = _modernize_array_constructor_delimiters(code)
    if "go" in lowered:
        code = GO_TO.sub("goto", code)

    # One combined substitution replaces the former loop of thirteen re.sub
    # calls per code segment.
    code = MULTIWORD_KEYWORD_SPACING.sub(lambda match: " ".join(match.group().lower().split()), code)
    code = COMPOUND_KEYWORD.sub(lambda match: match.group(1) + COMPOUND_KEYWORDS[match.group("keyword").lower()], code)

    if "end" in lowered:
        code = re.sub(r"(?<![A-Za-z0-9_])end[ \t]+", "end ", code, flags=re.IGNORECASE)
        code = re.sub(r"(?<![A-Za-z0-9_])end ([A-Za-z]+)[ \t]+", r"end \1 ", code, flags=re.IGNORECASE)
    if "do" in lowered:
        code = re.sub(r"(?<![A-Za-z0-9_])do[ \t]+", "do ", code, flags=re.IGNORECASE)
        code = re.sub(r"(?<![A-Za-z0-9_])do[ \t]+while[ \t]*\(", "do while (", code, flags=re.IGNORECASE)
    if "dimension" in lowered and "(" in code:
        code = re.sub(r"(?<![A-Za-z0-9_])dimension[ \t]*\(", "dimension(", code, flags=re.IGNORECASE)
    if "if" in lowered and "(" in code:
        code = re.sub(r"(?<![A-Za-z0-9_])if[ \t]*\(", "if (", code, flags=re.IGNORECASE)
    if "associate" in lowered and "(" in code:
        code = re.sub(r"(?<![A-Za-z0-9_])associate[ \t]*\(", "associate(", code, flags=re.IGNORECASE)
    if "result" in lowered and "(" in code:
        code = re.sub(r"(?<![A-Za-z0-9_])result[ \t]*\(", "result(", code, flags=re.IGNORECASE)
    if ("type" in lowered or "class" in lowered) and "(" in code:
        code = re.sub(
            r"(?<![A-Za-z0-9_])(?:type|class)[ \t]*\(",
            lambda match: match.group()[:-1].rstrip() + "(",
            code,
            flags=re.IGNORECASE,
        )
    if "select" in lowered and "type" in lowered and "(" in code:
        code = re.sub(r"\bselect[ \t]+type[ \t]*\(", "select type (", code, flags=re.IGNORECASE)
        code = re.sub(r"\bselect[ \t]+type[ \t]+is[ \t]*\(", "select type is (", code, flags=re.IGNORECASE)
    if "(" in code:
        code = PARENTHESIZED_STATEMENT.sub(lambda match: match.group()[:-1].rstrip() + "(", code)
    if "subroutine" in lowered and "(" in code:
        code = EMPTY_SUBROUTINE_ARGUMENTS.sub(r"\1", code)
    if "only" in lowered and ":" in code:
        code = re.sub(r"\bonly[ \t]*:", "only:", code, flags=re.IGNORECASE)
    if any(char in code for char in "([]"):
        code = re.sub(r"([([])[ \t]+", r"\1", code)
    if any(char in code for char in ")]"):
        code = re.sub(r"(?<=\S)[ \t]+([)\]])", r"\1", code)
    if ")" in code and "then" in lowered:
        code = re.sub(r"\)[ \t]*then\b", ") then", code, flags=re.IGNORECASE)
    if "(" in code and any(word in lowered for word in ("change", "form", "select", "sync")):
        code = re.sub(
            r"\b(?:change|form|select|sync)[ \t]+(?:rank|team)[ \t]*\(",
            lambda match: match.group()[:-1].rstrip() + " (",
            code,
            flags=re.IGNORECASE,
        )

    # Arithmetic IF and single-line IF statements need a space between the
    # closing condition and the following statement.
    match = IF_STATEMENT_START.match(code)
    if match is None:
        return code
    depth = 0
    for index in range(match.end() - 1, len(code)):
        if code[index] == "(":
            depth += 1
        elif code[index] == ")":
            depth -= 1
            if depth == 0:
                suffix = code[index + 1 :]
                if suffix:
                    stripped = suffix.lstrip()
                    if stripped and stripped != "&" and not stripped.lower().startswith("then"):
                        return code[: index + 1] + " " + stripped
                return code
    return code


def normalize_keyword_spacing(line: str, quote: str | None = None) -> tuple[str, str | None]:
    """Normalize compound keywords and spacing around Fortran syntax delimiters."""
    # The overwhelming common case contains no literal or comment, so avoid the
    # segment-building loop entirely.
    if quote is None and "'" not in line and '"' not in line and "!" not in line:
        return _normalize_keyword_spacing_code(line), None

    output: list[str] = []
    segment_start = 0
    index = 0
    while index < len(line):
        char = line[index]
        if quote:
            if char == quote:
                if line[index + 1 : index + 2] == quote:
                    index += 2
                    continue
                output.append(line[segment_start : index + 1])
                quote = None
                segment_start = index + 1
            index += 1
        elif char in "\"'":
            output.append(_normalize_keyword_spacing_code(line[segment_start:index]))
            quote = char
            segment_start = index
            index += 1
        elif char == "!":
            output.append(_normalize_keyword_spacing_code(line[segment_start:index]))
            output.append(line[index:])
            return "".join(output), quote
        else:
            index += 1
    output.append(line[segment_start:] if quote else _normalize_keyword_spacing_code(line[segment_start:]))
    return "".join(output), quote


def normalize_write_output_spacing(line: str, quote: str | None = None) -> str:
    """Put a space between a WRITE control list and an adjacent output item."""
    if "write" not in line.lower():
        return line
    insertions: list[int] = []
    index = 0
    active_quote = quote
    while index < len(line):
        char = line[index]
        if active_quote is not None:
            if char == active_quote:
                if line[index + 1 : index + 2] == active_quote:
                    index += 2
                    continue
                active_quote = None
            index += 1
            continue
        if char in "\"'":
            active_quote = char
            index += 1
            continue
        if char == "!":
            break
        if not (char.isalpha() or char == "_"):
            index += 1
            continue
        word_end = index + 1
        while word_end < len(line) and (line[word_end].isalnum() or line[word_end] == "_"):
            word_end += 1
        if line[index:word_end].lower() != "write":
            index = word_end
            continue

        opening = word_end
        while opening < len(line) and line[opening] in " \t":
            opening += 1
        if opening >= len(line) or line[opening] != "(":
            index = word_end
            continue

        depth = 1
        control_quote: str | None = None
        closing = opening + 1
        while closing < len(line) and depth:
            control_char = line[closing]
            if control_quote is not None:
                if control_char == control_quote:
                    if line[closing + 1 : closing + 2] == control_quote:
                        closing += 2
                        continue
                    control_quote = None
            elif control_char in "\"'":
                control_quote = control_char
            elif control_char == "(":
                depth += 1
            elif control_char == ")":
                depth -= 1
            closing += 1
        if depth:
            index = word_end
            continue

        if closing < len(line) and line[closing] not in " \t&!;\n":
            insertions.append(closing)
        index = closing

    if not insertions:
        return line
    insertion_points = set(insertions)
    return "".join(" " + char if index in insertion_points else char for index, char in enumerate(line))


_DECLARATION_START_WORDS = frozenset(
    {
        "integer",
        "real",
        "complex",
        "logical",
        "character",
        "type",
        "class",
        "procedure",
        "dimension",
        "allocatable",
        "pointer",
        "target",
        "optional",
        "parameter",
        "save",
        "value",
        "volatile",
        "asynchronous",
        "contiguous",
        "codimension",
    }
)


def _is_declaration_statement(code: str) -> bool:
    index = 0
    while index < len(code) and code[index].isspace():
        index += 1
    if index >= len(code) or not _is_identifier_start(code[index]):
        return False
    end = _identifier_end(code, index)
    word = code[index:end].lower()
    if word in _DECLARATION_START_WORDS:
        return True
    if word != "double":
        return False
    precision = end
    if precision >= len(code) or not code[precision].isspace():
        return False
    while precision < len(code) and code[precision].isspace():
        precision += 1
    if precision >= len(code) or not _is_identifier_start(code[precision]):
        return False
    precision_end = _identifier_end(code, precision)
    return code[precision:precision_end].lower() == "precision"


def normalize_delimiter_spacing(line: str, quote: str | None = None) -> tuple[str, str | None]:
    """Normalize spacing around declaration separators and non-quoted commas."""
    newline = "\n" if line.endswith("\n") else ""
    content = line.removesuffix(newline)
    comment_start = inline_comment_start(content, quote) if "!" in content else None
    code_end = comment_start if comment_start is not None else len(content)
    code = content[:code_end]

    if _is_declaration_statement(code):
        if "::" not in code:
            code = normalize_unseparated_declaration_spacing(code, quote)
        else:
            code = normalize_declaration_attribute_order(code)
        content = code + content[code_end:]
        code_end = len(code)

    # Most normalized lines contain neither delimiter. If there is no active or
    # newly opened literal, there is also no quote state to update.
    code = content[:code_end]
    if quote is None and "'" not in content and '"' not in content and "," not in code and "::" not in code:
        return content + newline, None

    output: list[str] = []
    index = 0
    while index < len(content):
        char = content[index]
        if quote:
            output.append(char)
            if char == quote:
                if index + 1 < len(content) and content[index + 1] == quote:
                    output.append(content[index + 1])
                    index += 2
                    continue
                quote = None
            index += 1
            continue

        if char in "\"'":
            quote = char
            output.append(char)
            index += 1
        elif char == "!":
            output.append(content[index:])
            break
        elif content.startswith("::", index):
            output.append("::")
            following = index + 2
            if following < len(content) and not content[following].isspace():
                output.append(" ")
            index += 2
        elif char == ",":
            indent = pop_trailing_whitespace(output)
            if not output:
                output.append(indent)
            output.append(",")
            following = index + 1
            while following < len(content) and content[following].isspace():
                following += 1
            if following < len(content):
                output.append(" ")
                index = following
            else:
                index += 1
        else:
            output.append(char)
            index += 1
    return "".join(output) + newline, quote


def remove_redundant_nested_parentheses(source: str) -> str:
    """Remove redundant nested ``((expression))`` parentheses.

    This is limited to ordinary right-hand-side expressions and ``if``/``while``
    conditions. Parentheses in procedure arguments and association contexts can
    affect whether an actual argument is definable, so they are left unchanged.
    Strings, comments, and preprocessor lines are ignored so source text is
    preserved.
    """

    def preceding_identifier(index: int, start: int = 0) -> str | None:
        """Return the identifier immediately before *index*, ignoring whitespace."""
        end = index
        while end > start and source[end - 1].isspace():
            end -= 1

        token_start = end
        while token_start > start and (source[token_start - 1].isalnum() or source[token_start - 1] == "_"):
            token_start -= 1

        if token_start == end:
            return None

        # Preserve the existing regex's behaviour for odd inputs such as
        # ``123name(``, where it finds the identifier suffix ``name``.
        match = re.search(rf"({IDENTIFIER})$", source[token_start:end])
        return None if match is None else match.group(1)

    removals: set[int] = set()

    # opening position, protected, safe, directly nested in parent
    stack: list[tuple[int, bool, bool, bool]] = []

    quote: str | None = None
    comment = False
    line_start = True
    physical_line_start = 0
    line_has_assignment = False

    # Last non-whitespace character. This lets us determine whether
    # parent + child openings are separated only by whitespace without
    # slicing and stripping source[outer + 1 : inner].
    last_nonwhitespace = -1

    index = 0
    length = len(source)

    while index < length:
        char = source[index]

        if comment:
            if char == "\n":
                comment = False
                line_start = True
                physical_line_start = index + 1
                line_has_assignment = False
            index += 1
            continue

        if line_start and char in " \t":
            index += 1
            continue

        if line_start and char == "#":
            comment = True
            line_start = False
            last_nonwhitespace = index
            index += 1
            continue

        if quote is not None:
            if char == quote:
                if source[index + 1 : index + 2] == quote:
                    index += 2
                    continue
                quote = None

            if char == "\n":
                physical_line_start = index + 1
                line_has_assignment = False

            index += 1
            line_start = False
            continue

        if char in "'\"":
            quote = char
            last_nonwhitespace = index

        elif char == "!":
            comment = True
            last_nonwhitespace = index

        else:
            # Equivalent to finding a non-"=>" ASSIGNMENT_BREAK_OPERATOR
            # in code_context(source[physical_line_start:index]).
            if char == "=":
                previous = source[index - 1] if index > physical_line_start else ""
                following = source[index + 1 : index + 2]

                if previous not in "<>=/" and following not in {"=", ">"}:
                    line_has_assignment = True

            if char == "(":
                preceding_name = preceding_identifier(index)
                condition_name = preceding_identifier(index, physical_line_start)

                is_argument_list = bool(preceding_name and preceding_name.lower() not in FORTRAN_STANDARD_WORDS)
                is_condition = bool(condition_name and condition_name.lower() in {"if", "while"})

                if stack:
                    parent_open, parent_protected, parent_safe, _ = stack[-1]
                    directly_nested = last_nonwhitespace == parent_open
                else:
                    parent_protected = False
                    parent_safe = False
                    directly_nested = False

                stack.append(
                    (
                        index,
                        is_argument_list or parent_protected,
                        line_has_assignment or is_condition or parent_safe,
                        directly_nested,
                    )
                )

            elif char == ")" and stack:
                inner, _, _, directly_nested = stack.pop()

                if stack:
                    _, protected, safe, _ = stack[-1]

                    if safe and not protected and directly_nested:
                        following = index + 1
                        while following < length and source[following].isspace():
                            following += 1

                        if following < length and source[following] == ")":
                            removals.update((inner, index))

            if not char.isspace():
                last_nonwhitespace = index

        index += 1

        if char == "\n":
            line_start = True
            physical_line_start = index
            line_has_assignment = False
        else:
            line_start = False

    if not removals:
        return source

    return "".join(char for index, char in enumerate(source) if index not in removals)


def normalize_unseparated_declaration_spacing(code: str, quote: str | None = None) -> str:
    """Collapse repeated spaces in an old-style declaration without touching literals."""
    code = re.sub(
        r"^(\s*(?:integer|real|complex|logical|character|type|class)\s*\([^)]*\))(?=[A-Za-z_])",
        r"\1 ",
        code,
        count=1,
        flags=re.IGNORECASE,
    )
    leading_length = len(code) - len(code.lstrip(" \t"))
    output: list[str] = [code[:leading_length]]
    pending_space = False
    for char in code[leading_length:]:
        if quote:
            output.append(char)
            if char == quote:
                quote = None
            continue
        if char in "\"'":
            if pending_space:
                output.append(" ")
                pending_space = False
            quote = char
            output.append(char)
        elif char in " \t":
            pending_space = True
        else:
            if pending_space:
                output.append(" ")
                pending_space = False
            output.append(char)
    return "".join(output)


def normalize_declaration_attribute_order(code: str) -> str:
    """Place an OPTIONAL declaration attribute after the other attributes."""
    separator = code.find("::")
    if separator < 0:
        return code
    specification = code[:separator]
    attributes = _split_top_level(specification)
    optional = [attribute for attribute in attributes if attribute.strip().lower() == "optional"]
    if not optional:
        return code
    attributes = [attribute for attribute in attributes if attribute.strip().lower() != "optional"]
    return ",".join(attributes + optional) + code[separator:]


def declaration_separator_info(line: str) -> tuple[int, int, int] | None:
    """Return the column and adjacent whitespace counts for a code ``::``."""
    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote:
            if char == quote:
                if line[index + 1 : index + 2] == quote:
                    index += 2
                    continue
                quote = None
            index += 1
        elif char in "\"'":
            quote = char
            index += 1
        elif char == "!":
            return None
        elif line.startswith("::", index):
            before = index
            while before > 0 and line[before - 1] in " \t":
                before -= 1
            after = index + 2
            while after < len(line) and line[after] in " \t":
                after += 1
            return index, index - before, after - (index + 2)
        else:
            index += 1
    return None


def normalize_declaration_separator_alignment(lines: list[str]) -> list[str]:
    """Compress adjacent declaration alignment without ever adding padding."""
    cpp_lines: list[bool] = []
    cpp_continuation = False
    for line in lines:
        cpp_line = cpp_continuation or is_preprocessor_line(line)
        cpp_lines.append(cpp_line)
        cpp_continuation = cpp_line and cpp_line_continues(line)
    separators = [None if cpp_line else declaration_separator_info(line) for line, cpp_line in zip(lines, cpp_lines)]
    updated = lines.copy()

    def normalize_block(block: list[tuple[int, tuple[int, int, int]]]) -> None:
        if not block:
            return
        separator_column = max(column - before + 1 for _, (column, before, _) in block)
        can_compress_alignment = len(block) > 1 and all(column >= separator_column for _, (column, _, _) in block)
        for line_index, (column, before, after) in block:
            prefix_end = column - before
            target_column = separator_column if can_compress_alignment else prefix_end + 1
            suffix_start = column + 2 + after
            before_spaces = " " * (target_column - prefix_end)
            updated[line_index] = (
                lines[line_index][:prefix_end] + before_spaces + ":: " + lines[line_index][suffix_start:]
            )

    index = 0
    while index < len(lines):
        if separators[index] is None:
            index += 1
            continue

        block_indices: list[int] = []
        while index < len(lines):
            if separators[index] is not None:
                block_indices.append(index)
                index += 1
                continue
            if lines[index].lstrip().startswith("!") and not cpp_lines[index]:
                index += 1
                continue
            break
        if not block_indices:
            index += 1
            continue
        block: list[tuple[int, tuple[int, int, int]]] = []
        separator_column = 0
        for line_index in block_indices:
            separator = separators[line_index]
            minimum_column = separator[0] - separator[1] + 1
            proposed_column = max(separator_column, minimum_column)
            if block and (
                separator[0] < proposed_column or any(column < proposed_column for _, (column, _, _) in block)
            ):
                normalize_block(block)
                block = []
                separator_column = 0
                proposed_column = minimum_column
            block.append((line_index, separator))
            separator_column = proposed_column
        normalize_block(block)
    return updated


def limit_blank_lines(lines: list[str], maximum: int = 2) -> list[str]:
    """Limit each run of whitespace-only lines to *maximum* lines."""
    limited: list[str] = []
    blank_count = 0
    cpp_continuation = False
    for line in lines:
        cpp_line = cpp_continuation or is_preprocessor_line(line)
        if cpp_line:
            limited.append(line)
            cpp_continuation = cpp_line_continues(line)
            blank_count = 0
            continue
        if line.strip():
            blank_count = 0
            limited.append(line)
        elif blank_count < maximum:
            blank_count += 1
            limited.append(line)
    return limited


def normalize_program_unit_spacing(lines: list[str]) -> list[str]:
    """Leave one blank line around program-unit boundaries without touching CPP bodies."""
    normalized: list[str] = []
    unit_depth = 0
    type_depth = 0
    interface_depth = 0
    add_blank_before_next = False
    states = scan_physical_lines(lines)
    for state in states:
        line = state.text
        if state.is_cpp:
            if add_blank_before_next and not state.cpp_continuation_in:
                if normalized and normalized[-1].strip():
                    normalized.append("\n")
                add_blank_before_next = False
            normalized.append(line)
            continue

        code = code_context(line)
        is_blank = state.is_blank
        if add_blank_before_next:
            if is_blank:
                continue
            normalized.append("\n")
            add_blank_before_next = False

        if is_blank and (unit_depth > 0 or interface_depth > 0) and normalized and not normalized[-1].strip():
            continue

        if TYPE_DEFINITION_END.match(code):
            type_depth = max(0, type_depth - 1)
        elif TYPE_DEFINITION_START.match(code):
            type_depth += 1

        is_end = interface_depth == 0 and (
            PROGRAM_UNIT_END.match(code) is not None or BARE_PROGRAM_UNIT_END.match(code) is not None
        )
        header = None if is_end else _scope_header_names(code)
        if interface_depth == 0 and (header or MODULE_DECLARATION.match(code)):
            unit_depth += 1

        is_contains = unit_depth > 0 and type_depth == 0 and CONTAINS_STATEMENT.match(code)
        if is_contains or is_end:
            while normalized and not normalized[-1].strip():
                normalized.pop()
            if normalized:
                normalized.append("\n")

        normalized.append(line)
        if is_contains:
            add_blank_before_next = True
        if is_end:
            unit_depth = max(0, unit_depth - 1)
        if INTERFACE_END.match(code):
            interface_depth = max(0, interface_depth - 1)
        elif INTERFACE_START.match(code):
            interface_depth += 1
    return normalized


def normalize_output_whitespace(text: str) -> str:
    """Remove trailing spaces and ensure one final newline, like pre-commit."""
    if not text:
        return text
    lines: list[str] = []
    cpp_continuation = False
    for line in text.split("\n"):
        cpp_line = cpp_continuation or is_preprocessor_line(line)
        lines.append(line if cpp_line else line.rstrip(" \t"))
        cpp_continuation = cpp_line and cpp_line_continues(line)
    while len(lines) > 1 and not lines[-1]:
        lines.pop()
    return "\n".join(lines) + "\n"


def dominant_line_ending(text: str) -> str:
    """Return the most common newline sequence in *text*."""
    endings = re.findall(r"\r\n|\r|\n", text)
    return Counter(endings).most_common(1)[0][0] if endings else "\n"


def comment_block_lines(lines: list[str]) -> set[int]:
    """Return full-line comments in consecutive blocks where alignment matters."""

    def is_full_comment(line: str) -> bool:
        content = line.lstrip(" \t")
        return content.startswith("!") and not content.startswith("!!") and DIRECTIVE_SENTINEL.match(content) is None

    preserved: set[int] = set()
    index = 0
    while index < len(lines):
        if not is_full_comment(lines[index]):
            index += 1
            continue
        block_start = index
        while index < len(lines) and is_full_comment(lines[index]):
            index += 1
        if index - block_start > 1:
            preserved.update(range(block_start, index))
    return preserved


def normalize_comment_spacing(line: str, quote: str | None = None, preserve_after: bool = False) -> str:
    """Ensure ordinary comments have one space after ``!`` and before inline ``!``."""
    comment_start = inline_comment_start(line, quote)
    if (
        comment_start is None
        or line[comment_start : comment_start + 2] == "!!"
        or DIRECTIVE_SENTINEL.match(line, comment_start)
    ):
        return line

    before = line[:comment_start]
    leading = before[: len(before) - len(before.lstrip())]
    code = before[len(leading) :].rstrip()
    comment = line[comment_start + 1 :]
    if not comment.strip():
        return line
    prefix = leading + (code + " " if code else "")
    if preserve_after and comment[:1] in " \t":
        return prefix + "!" + comment
    return prefix + "! " + comment.lstrip(" \t")


def comment_contains_assignment(comment: str) -> bool:
    """Return whether a comment starts with a plausible Fortran assignment."""
    return COMMENTED_ASSIGNMENT.match(comment) is not None


def format_comment_operators(comment: str) -> str:
    """Normalize operators in comments while preserving quoted examples."""
    output: list[str] = []
    quote: str | None = None
    index = 0
    while index < len(comment):
        char = comment[index]
        if quote:
            output.append(char)
            if char == quote:
                quote = None
            index += 1
            continue
        if char in "\"'":
            quote = char
            output.append(char)
            index += 1
            continue
        operator = LEGACY_OPERATOR.match(comment, index)
        if operator:
            append_spaced_operator(output, MODERN_OPERATOR[operator.group(1).lower()])
            index = skip_operator_whitespace(comment, operator.end())
            continue
        operator = SPACED_OPERATOR.match(comment, index)
        if operator:
            if operator.group() == "=" and is_named_parameter(comment, index):
                append_compact_operator(output, operator.group())
            else:
                append_spaced_operator(output, operator.group())
            index = skip_operator_whitespace(comment, operator.end())
            continue
        operator = ARITHMETIC_OPERATOR.match(comment, index)
        if operator and operator.group() == "+" and is_binary_arithmetic_operator(comment, index, operator.group()):
            append_spaced_operator(output, operator.group())
            index = skip_operator_whitespace(comment, operator.end())
            continue
        output.append(char)
        index += 1
    return "".join(output).replace(" \n", "\n")


def inline_comment_start(line: str, quote: str | None = None) -> int | None:
    """Return the position of an unquoted Fortran comment marker.

    *quote* is the character-literal delimiter left open by the preceding line.
    """
    if quote is None and "'" not in line and '"' not in line:
        index = line.find("!")
        return None if index < 0 else index
    for index, char in enumerate(line):
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "!":
            return index
    return None


def quote_after(line: str, quote: str | None = None) -> str | None:
    """Return the character-literal delimiter *line* leaves open, if any.

    *quote* is the delimiter left open by the preceding line. A literal that is
    still open at the end of a line continues into the next one, which makes
    that line's leading `&` and the spacing before its trailing `&` part of the
    literal rather than layout that may be rewritten.
    """
    content = line.rstrip("\n")
    if quote is None and "'" not in content and '"' not in content:
        return None
    index = 0
    while index < len(content):
        char = content[index]
        if quote:
            if char == quote:
                if content[index + 1 : index + 2] == quote:
                    index += 2
                    continue
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "!":
            break
        index += 1
    return quote


def blank_leading_continuation(line: str) -> str:
    """Return *line* with a leading continuation marker replaced by a space.

    The marker is optional in free-form Fortran and carries no meaning outside
    a character literal, so blanking it keeps every other column in place.
    This is used while matching identifiers, where preserving source columns
    is necessary for context offsets.
    """
    content = line.rstrip("\n")
    newline = line[len(content) :]
    indent_length = len(content) - len(content.lstrip())
    if content[indent_length : indent_length + 1] != "&":
        return line
    return content[:indent_length] + " " + content[indent_length + 1 :] + newline


def remove_leading_continuation(line: str) -> str:
    """Remove a leading continuation marker and following layout whitespace."""
    content = line.rstrip("\n")
    newline = line[len(content) :]
    indent_length = len(content) - len(content.lstrip())
    if content[indent_length : indent_length + 1] != "&":
        return line
    return content[:indent_length] + content[indent_length + 1 :].lstrip(" \t") + newline


def normalize_continuation_marker(line: str) -> str:
    """Ensure an ending free-form continuation marker has one preceding space."""
    content = line.rstrip("\n")
    newline = line[len(content) :]
    comment_start = inline_comment_start(content)
    code_end = comment_start if comment_start is not None else len(content)
    code = content[:code_end]
    if not code.rstrip().endswith("&"):
        return line
    stripped_code = code.rstrip()
    suffix = code[len(stripped_code) :] if comment_start is not None else ""
    normalized = stripped_code[:-1].rstrip() + " &" + suffix
    return normalized + content[code_end:] + newline


def mask_quoted_text(text: str) -> str:
    """Replace quoted spans with filler so operator scans never match inside them."""
    if "'" not in text and '"' not in text:
        return text
    output: list[str] = []
    quote: str | None = None
    for char in text:
        if quote:
            output.append("#")
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
            output.append("#")
        else:
            output.append(char)
    return "".join(output)


def paren_depth_profile(text: str) -> list[int]:
    """Return the bracket-nesting depth before each index of *text*, outside quotes."""
    depths = [0]
    depth = 0
    quote: str | None = None
    for char in text:
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
        depths.append(depth)
    return depths


def scan_limit(line: str, limit: int) -> int:
    """Restrict *limit* to stop before an unquoted comment marker, if any."""
    comment_start = inline_comment_start(line)
    return min(limit, comment_start) if comment_start is not None else limit


def operator_break_position(masked: str, depths: list[int], start: int, limit: int) -> int | None:
    """Return the best tiered break in ``masked[start:limit]``, if there is one.

    The search runs over the whole line and filters by position afterwards:
    passing *limit* to finditer() as an end position would hide the character
    after it from the patterns' lookaheads, so a `/=` or `==` straddling the
    limit would look like a lone `/` or `=` and be split down the middle.
    """
    minimum_fill = int(limit * MINIMUM_BREAK_FILL)
    best: tuple[bool, int, int, int] | None = None
    for tier, pattern in enumerate(BREAK_OPERATORS):
        for match in pattern.finditer(masked, start):
            position = match.end()
            if position > limit:
                break
            candidate = (position < minimum_fill, depths[position], tier, -position)
            if best is None or candidate < best:
                best = candidate
    return None if best is None else -best[3]


def assignment_break_position(line: str, masked: str, limit: int) -> int | None:
    """Return the last top-level ``=``/``=>`` break before *limit*, if there is one."""
    position: int | None = None
    for match in ASSIGNMENT_BREAK_OPERATOR.finditer(masked):
        if match.end() > limit:
            break
        if match.group() == "=" and is_named_parameter(line, match.start()):
            continue
        position = match.end()
    return position


def statement_head_end(line: str, masked: str, depths: list[int], limit: int) -> int:
    """Return the end of a statement's head, i.e. what a break should not split.

    The head is a declaration's attribute list up to ``::``, or an assignment's
    left-hand side up to ``=``. Breaking inside either splits a name or an
    attribute list off from what it belongs to, when the boundary right after
    the head is both available and the one a human would use.
    """
    declaration = masked.find("::", 0, limit)
    if declaration >= 0 and depths[declaration] == 0:
        return declaration + 2
    for match in ASSIGNMENT_BREAK_OPERATOR.finditer(masked):
        if match.end() > limit:
            break
        if depths[match.start()] == 0 and not is_named_parameter(line, match.start()):
            return match.end()
    return 0


def wrap_position(line: str, limit: int) -> int | None:
    """Return a safe wrapping position before *limit* outside quoted text.

    Candidates are ranked by the shallowest bracket-nesting depth first, then
    by how loosely the operator binds (loosest first), so a break lands on the
    same outermost boundary a human would choose rather than the last operator
    that happens to fit before the column limit.

    A statement's head — a declaration's attribute list up to ``::``, or an
    assignment's left-hand side up to ``=`` — is left intact whenever the
    boundary at its end is within reach, since that boundary is the one a
    human would use: breaks are looked for after it, and it is taken as the
    break when there are none. Only a head that overruns *limit* by itself is
    searched inside.
    """
    limit = scan_limit(line, limit)
    masked = mask_quoted_text(line)
    depths = paren_depth_profile(line)

    head_end = statement_head_end(line, masked, depths, limit)
    if head_end:
        position = operator_break_position(masked, depths, head_end, limit)
        if position is not None:
            return position
        return head_end

    position = operator_break_position(masked, depths, 0, limit)
    if position is not None:
        return position

    position = assignment_break_position(line, masked, limit)
    if position is not None:
        return position

    whitespace_position: int | None = None
    quote = None
    for index, char in enumerate(line[:limit]):
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char.isspace():
            whitespace_position = index
    return whitespace_position


def assignment_wrap_position(line: str, limit: int) -> int | None:
    """Return an assignment break point that is not a comparison or named argument."""
    quote: str | None = None
    for index, char in enumerate(line[:limit]):
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "!":
            return None
        elif (
            char == "="
            and (index == 0 or line[index - 1] not in "<>/=")
            and (index + 1 == len(line) or line[index + 1] not in "=>")
            and not is_named_parameter(line, index)
        ):
            return index + 1
    return None


def wrap_line(line: str, line_length: int = MAX_LINE_LENGTH) -> list[str]:
    """Wrap one free-form Fortran code line using continuation markers."""
    newline = "\n" if line.endswith("\n") else ""
    content = line.removesuffix(newline)
    if len(content) <= line_length or is_preprocessor_line(line) or content.lstrip().startswith("!"):
        return [line]

    indent_length = len(content) - len(content.lstrip())
    indent = content[:indent_length]
    comment_start = inline_comment_start(content)
    if comment_start is not None:
        # A comment cannot sit between a continuation marker and the text it
        # continues, so lift it above the statement before wrapping the code.
        code = content[:comment_start].rstrip()
        return [indent + content[comment_start:].lstrip() + newline, *wrap_line(code + newline, line_length)]

    body = content[indent_length:].rstrip()
    continuation = body.endswith("&")
    if continuation:
        body = body[:-1].rstrip()

    wrapped: list[str] = []
    suffix = " &" if continuation else ""
    continuation_indent = indent + "    "
    current_indent = indent
    first_break = True
    while len(current_indent) + len(body) + len(suffix) > line_length:
        position = None
        if first_break:
            assignment_position = assignment_wrap_position(body, line_length - len(current_indent) - len(" &"))
            if assignment_position is not None:
                remainder = body[assignment_position:].lstrip()
                if len(continuation_indent) + len(remainder) + len(suffix) <= line_length:
                    position = assignment_position
        if position is None:
            position = wrap_position(body, line_length - len(current_indent) - len(" &"))
        if position is None:
            return [line]
        text = (current_indent + body[:position]).rstrip()
        wrapped.append(text + " &\n")
        body = body[position:].lstrip()
        current_indent = continuation_indent
        first_break = False
    wrapped.append(current_indent + body + suffix + newline)
    return wrapped


OPENMP_DIRECTIVE = re.compile(r"^(?P<indent>[ \t]*)!\$omp(?P<continuation>&?)(?P<body>.*)$", re.IGNORECASE)


def join_openmp_directive(lines: list[str]) -> str | None:
    """Join a free-form OpenMP directive and its continuation lines."""
    parts: list[str] = []
    indent = ""
    for index, line in enumerate(lines):
        match = OPENMP_DIRECTIVE.match(line.rstrip("\n"))
        if match is None or (not index and match.group("continuation")):
            return None
        if index == 0:
            indent = match.group("indent")
        body = match.group("body").lstrip()
        if index and body.startswith("&"):
            body = body[1:].lstrip()
        if inline_comment_start(body) is not None:
            return None
        if index < len(lines) - 1:
            if not body.endswith("&"):
                return None
            body = body[:-1].rstrip()
        parts.append(body)
    newline = "\n" if lines[-1].endswith("\n") else ""
    return indent + "!$OMP " + " ".join(parts) + newline


def wrap_openmp_directive(line: str, line_length: int = MAX_LINE_LENGTH) -> list[str]:
    """Wrap a free-form OpenMP directive with repeated ``!$OMP`` sentinels."""
    newline = "\n" if line.endswith("\n") else ""
    content = line.removesuffix(newline)
    match = OPENMP_DIRECTIVE.match(content)
    if match is None or match.group("continuation") or len(content) <= line_length:
        return [line]

    indent = match.group("indent")
    body = match.group("body").strip()
    if inline_comment_start(body) is not None:
        return [line]

    wrapped: list[str] = []
    prefix = indent + "!$OMP "
    continuation_prefix = indent + "!$OMP "
    while len(prefix) + len(body) > line_length:
        position = wrap_position(body, line_length - len(prefix) - len(" &"))
        if position is None:
            return [line]
        wrapped.append(prefix + body[:position].rstrip() + " &\n")
        body = body[position:].lstrip()
        prefix = continuation_prefix
    wrapped.append(prefix + body + newline)
    return wrapped


def ends_with_continuation(line: str, quote: str | None = None) -> bool:
    """Return whether a line's code portion ends with a continuation marker."""
    content = line.rstrip("\n")
    openmp = OPENMP_DIRECTIVE.match(content)
    if openmp is not None:
        return openmp.group("body").rstrip().endswith("&")
    comment_start = inline_comment_start(content, quote)
    if comment_start is not None:
        content = content[:comment_start]
    return content.rstrip().endswith("&")


def statement_context(prefix: str, line: str, quote: str | None = None) -> str:
    """Return a statement's code through *line*, as context for its next line.

    Comments and continuation markers are dropped, and a space separates the
    lines so that the last token of one never runs into the first of the next.
    """
    content = line.rstrip("\n")
    comment_start = inline_comment_start(content, quote)
    if comment_start is not None:
        content = content[:comment_start]
    content = content.rstrip().removesuffix("&").strip()
    if prefix:
        content = content.removeprefix("&")
    return prefix + content + " "


def continues_statement(line: str, quote: str | None = None) -> tuple[bool, str | None]:
    """Return whether *line* is continued, and the character literal it leaves open.

    *quote* is the literal delimiter left open by the preceding line. An open
    literal always continues, since a literal cannot end at a line break.
    """
    if quote is None and "'" not in line and '"' not in line:
        content = line.rstrip("\n")
        openmp = OPENMP_DIRECTIVE.match(content)
        if openmp is not None:
            return openmp.group("body").rstrip().endswith("&"), None
        comment_start = content.find("!")
        if comment_start >= 0:
            content = content[:comment_start]
        return content.rstrip().endswith("&"), None
    closing_quote = quote_after(line, quote)
    return closing_quote is not None or ends_with_continuation(line, quote), closing_quote


def join_continued_lines(lines: list[str]) -> str:
    """Join a continued Fortran statement into one physical line."""
    first = lines[0].rstrip("\n")
    indent_length = len(first) - len(first.lstrip())
    indent = first[:indent_length]
    parts: list[str] = []

    for index, line in enumerate(lines):
        content = line.rstrip("\n").lstrip()
        if index > 0 and content.startswith("&"):
            content = content[1:].lstrip()
        if index < len(lines) - 1:
            content = content.rstrip()
            if not content.endswith("&"):
                raise ValueError("Expected a Fortran continuation marker")
            content = content[:-1].rstrip()
        parts.append(content)

    joined = indent + parts[0]
    for part in parts[1:]:
        separator = "" if joined.rstrip().endswith(("*", "/")) else " "
        joined += separator + part

    newline = "\n" if lines[-1].endswith("\n") else ""
    return joined + newline


def _lexical_continuation_prefix(line: str, quote: str | None = None) -> str | None:
    """Return the code before a trailing marker when it ends a token prefix."""
    if quote is not None or is_preprocessor_line(line) or quote_after(line, quote) is not None:
        return None
    content = line.rstrip("\n")
    comment_start = inline_comment_start(content, quote)
    if comment_start is not None:
        return None
    code = content.rstrip()
    if not code.endswith("&"):
        return None

    prefix = code[:-1]
    if not prefix or not (prefix[-1].isalnum() or prefix[-1] == "_"):
        return None
    return prefix


def _leading_continuation_suffix(line: str) -> str | None:
    """Return the text after a leading free-form continuation marker."""
    content = line.rstrip("\n")
    indent_length = len(content) - len(content.lstrip())
    if content[indent_length : indent_length + 1] != "&":
        return None
    suffix = content[indent_length + 1 :]
    if not suffix or not (suffix[0].isalnum() or suffix[0] == "_"):
        return None
    return suffix


def join_lexical_token_continuations(source: str) -> str:
    """Join free-form continuations that split a token before formatting it."""
    lines = source.splitlines(keepends=True)
    joined: list[str] = []
    index = 0
    quote: str | None = None
    cpp_continuation = False
    while index < len(lines):
        line = lines[index]
        if cpp_continuation or is_preprocessor_line(line):
            joined.append(line)
            cpp_continuation = cpp_line_continues(line)
            quote = None
            index += 1
            continue
        while index + 1 < len(lines):
            prefix = _lexical_continuation_prefix(line, quote)
            if prefix is None:
                break
            next_line = lines[index + 1]
            suffix = _leading_continuation_suffix(next_line)
            if suffix is None:
                break
            newline = "\n" if next_line.endswith("\n") else ""
            line = prefix + suffix + newline
            index += 1
        joined.append(line)
        quote = quote_after(line, quote)
        index += 1
    return "".join(joined)


def is_lexical_token_continuation(first_line: str, second_line: str) -> bool:
    """Return whether a continued boundary splits one lexical token."""
    return _trailing_token_continuation(first_line) and _leading_continuation_suffix(second_line) is not None


def normalize_continuations(lines: list[str]) -> list[str]:
    """Normalize trailing markers and replace leading continuation markers with spaces.

    Lines that resume a character literal are left untouched: their leading `&`
    is required to keep the literal's text intact, and the spaces before their
    trailing `&` belong to the literal.
    """
    normalized: list[str] = []
    continuation = False
    cpp_continuation = False
    quote: str | None = None
    for index, line in enumerate(lines):
        if cpp_continuation or is_preprocessor_line(line):
            normalized.append(line)
            cpp_continuation = cpp_line_continues(line)
            continuation = False
            quote = None
            continue

        opening_quote = quote
        continues, quote = continues_statement(line, opening_quote)
        in_literal = opening_quote is not None or quote is not None
        lexical_prefix = index > 0 and is_lexical_token_continuation(lines[index - 1], line)
        lexical_suffix = index + 1 < len(lines) and is_lexical_token_continuation(line, lines[index + 1])
        if continuation and not in_literal and not lexical_prefix:
            line = remove_leading_continuation(line)
        if not in_literal and not lexical_suffix:
            line = normalize_continuation_marker(line)
        normalized.append(line)
        continuation = continues
    return normalized


def normalize_openmp_continuation_sentinels(lines: list[str]) -> list[str]:
    """Use ``!$OMP`` for every line of a continued OpenMP directive."""
    normalized: list[str] = []
    continuation = False
    for line in lines:
        newline = "\n" if line.endswith("\n") else ""
        content = line.removesuffix(newline)
        match = OPENMP_DIRECTIVE.match(content)
        if match is None:
            normalized.append(line)
            continuation = False
            continue

        body = match.group("body")
        if match.group("continuation") or (continuation and body.lstrip().startswith("&")):
            body = body.lstrip().removeprefix("&").lstrip()
            line = match.group("indent") + "!$OMP " + body + newline
        normalized.append(line)
        continuation = ends_with_continuation(line)
    return normalized


def remove_terminal_procedure_returns(
    lines: list[str], procedure_cases: Iterable[ProcedureDeclarationCases] = ()
) -> list[str]:
    """Remove a bare final ``return`` from function and subroutine bodies."""
    source = "".join(lines)
    statements = _code_statements(source)
    statements_by_end_line = {statement.end_line: statement for statement in statements}
    statement_index_by_end_line = {statement.end_line: index for index, statement in enumerate(statements)}
    procedures = tuple(procedure_cases) or extract_procedure_cases(source)
    removed_lines: set[int] = set()
    for procedure in procedures:
        header = statements_by_end_line.get(procedure.start_line)
        end = statements_by_end_line.get(procedure.end_line)
        if header is None or end is None:
            continue
        header_match = FUNCTION_OR_SUBROUTINE_WORD.search(header.text)
        if header_match is None:
            continue
        end_index = statement_index_by_end_line.get(procedure.end_line)
        if end_index is None or end_index == 0:
            continue
        previous = statements[end_index - 1]
        if not (procedure.start_line < previous.end_line < procedure.end_line):
            continue
        if previous.text.lower() != "return" or previous.start_line != previous.end_line:
            continue
        line = lines[previous.start_line]
        if code_context(line).strip().lower() == "return" and inline_comment_start(line) is None:
            removed_lines.add(previous.start_line)
    if not removed_lines:
        return lines
    return [line for line_number, line in enumerate(lines) if line_number not in removed_lines]


def has_continued_string(lines: list[str]) -> bool:
    """Return whether a character literal crosses a physical source-line boundary."""
    quote: str | None = None
    for line in lines:
        quote = quote_after(line, quote)
        if quote is not None:
            return True
    return False


def detach_final_inline_comment(statement: list[str]) -> tuple[list[str], list[str]] | None:
    """Move a final inline comment above a statement before it is reflowed."""
    comments: list[tuple[int, int]] = []
    for index, line in enumerate(statement):
        comment_start = inline_comment_start(line.rstrip("\n"))
        if comment_start is not None:
            comments.append((index, comment_start))
    if not comments:
        return [], statement
    if len(comments) != 1 or comments[0][0] != len(statement) - 1:
        return None

    line_index, comment_start = comments[0]
    line = statement[line_index]
    newline = "\n" if line.endswith("\n") else ""
    comment_indent = line[: len(line) - len(line.lstrip())]
    comment = comment_indent + line[comment_start:].lstrip()
    if newline and not comment.endswith("\n"):
        comment += newline
    code_line = line[:comment_start].rstrip() + newline
    return [comment], [*statement[:-1], code_line]


def rewrap_lines(
    lines: list[str],
    line_length: int = MAX_LINE_LENGTH,
    procedure_cases: Iterable[ProcedureDeclarationCases] = (),
    preserved_names: Collection[str] = frozenset(),
    preprocessor_names: Collection[str] = frozenset(),
    scoped_declared_names: Iterable[ScopedDeclaredNames] = (),
    uppercase_single_l: bool = False,
) -> list[str]:
    """Reflow long statements using the shared physical-line scanner."""
    procedure_cases = tuple(procedure_cases)
    states = scan_physical_lines(lines)
    wrapped: list[str] = []
    index = 0
    while index < len(lines):
        state = states[index]
        if state.is_cpp or state.is_blank:
            wrapped.append(lines[index])
            index += 1
            continue
        if state.is_comment:
            openmp_directive = join_openmp_directive([lines[index]])
            if openmp_directive is not None and len(lines[index].rstrip("\n")) > line_length:
                wrapped.extend(wrap_openmp_directive(openmp_directive, line_length))
            else:
                wrapped.append(lines[index])
            index += 1
            continue

        end_index = _statement_end_index(states, index)
        statement = lines[index : end_index + 1]
        start_index = index
        index = end_index + 1

        openmp_directive = join_openmp_directive(statement)
        if openmp_directive is not None and any(len(line.rstrip("\n")) > line_length for line in statement):
            wrapped.extend(wrap_openmp_directive(openmp_directive, line_length))
            continue

        has_long_line = any(len(line.rstrip("\n")) > line_length for line in statement)
        detached = detach_final_inline_comment(statement)
        has_unsafe_content = any(
            states[pos].is_cpp or states[pos].is_comment or states[pos].is_blank
            for pos in range(start_index, end_index + 1)
        ) or has_continued_string(statement)
        if has_long_line and detached is not None and not has_unsafe_content:
            comments, code_statement = detached
            active_procedure = active_procedure_at(procedure_cases, start_index)
            local_names = active_procedure.local_names if active_procedure else frozenset()
            file_declared_names = declared_names_at(scoped_declared_names, start_index)
            joined, _ = lowercase_line(
                join_continued_lines(code_statement),
                preserved_names=preserved_names,
                local_names=local_names,
                preprocessor_names=preprocessor_names,
                file_declared_names=file_declared_names,
                uppercase_single_l=uppercase_single_l,
            )
            # Same order as the per-line pass in format_text(): lowercase_line
            # pads `.not.` as an operator, and only keyword spacing removes the
            # space it leaves after an opening bracket (`if ( .not. (`).
            joined, _ = normalize_keyword_spacing(joined)
            joined, _ = normalize_delimiter_spacing(joined)
            wrapped.extend(comments)
            wrapped.extend(wrap_line(joined, line_length))
        else:
            wrapped.extend(statement)
    return wrapped


def format_text(
    original: str,
    wrap: bool = True,
    module_cases: Mapping[str, str] | None = None,
    symbol_cases: Mapping[str, str] | None = None,
    procedure_cases: Iterable[ProcedureDeclarationCases] = (),
    scope_cases: Iterable[NamedScopeCase] = (),
    type_procedure_cases: Mapping[str, str] | None = None,
    type_component_cases: Mapping[tuple[str, str], str] | None = None,
    variable_type_cases: Mapping[str, str] | None = None,
    type_component_type_cases: Mapping[tuple[str, str], str] | None = None,
    uppercase_single_l: bool = False,
    macro_cases: Mapping[str, str] | Collection[str] = (),
) -> str:
    """Format Fortran source text without reading from or writing to a file."""
    line_ending = dominant_line_ending(original)
    source = original.replace("\r\n", "\n").replace("\r", "\n")
    explicit_macro_cases = normalize_macro_cases(macro_cases)
    source_macro_cases = extract_preprocessor_cases(source)
    # Source-defined macros keep their existing exact-match behavior. Only
    # explicit -D/--define inputs are authoritative enough to canonicalize case.
    source = replace_preprocessor_cases(source, explicit_macro_cases)
    procedure_cases = tuple(procedure_cases) or extract_procedure_cases(source)
    symbol_cases = symbol_cases or {}
    preprocessor_names = frozenset(source_macro_cases.values()) | frozenset(explicit_macro_cases.values())
    source = replace_declared_cases(
        source,
        module_cases or {},
        symbol_cases,
        procedure_cases,
        scope_cases,
        type_procedure_cases,
        type_component_cases,
        variable_type_cases,
        type_component_type_cases,
        preprocessor_names,
    )
    joined_source = join_lexical_token_continuations(source)
    lexical_lines_changed = joined_source != source
    source = remove_redundant_nested_parentheses(joined_source)
    # Only lexical token joins can collapse physical lines here. Casing and
    # parenthesis cleanup preserve line ranges and declaration spellings, so
    # retain the existing procedure analysis unless a join actually occurred.
    if lexical_lines_changed:
        procedure_cases = extract_procedure_cases(source)
    scoped_declared_names = extract_scoped_declared_names(source)
    prefix = ""
    lines: list[str] = []
    source_lines, source_states = _scan_source(source)
    active_procedures = _active_procedures_by_line(procedure_cases, len(source_states))
    declared_names_per_line = _declared_names_by_line(scoped_declared_names, len(source_states))
    preserved_comment_lines = comment_block_lines(source_lines)
    for state in source_states:
        line_number = state.number
        line = state.text
        if state.is_cpp:
            lines.append(line)
            prefix = ""
            continue
        if state.is_blank:
            lines.append(line)
            continue
        if state.is_comment and state.quote_in is not None:
            lines.append(line)
            continue

        starting_quote = state.quote_in
        active_procedure = active_procedures[line_number]
        local_names = active_procedure.local_names if active_procedure else frozenset()
        lowercase, _ = lowercase_line(
            line,
            starting_quote,
            prefix,
            preserved_names=symbol_cases,
            local_names=local_names,
            preprocessor_names=preprocessor_names,
            file_declared_names=declared_names_per_line[line_number],
            uppercase_single_l=uppercase_single_l,
        )
        normalized, _ = normalize_keyword_spacing(lowercase, starting_quote)
        normalized = normalize_write_output_spacing(normalized, starting_quote)
        if not is_preprocessor_line(line):
            normalized, _ = normalize_delimiter_spacing(normalized, starting_quote)
            normalized = normalize_comment_spacing(
                normalized,
                starting_quote,
                preserve_after=line_number in preserved_comment_lines,
            )
        lines.append(normalized)
        if not state.is_comment:
            prefix = statement_context(prefix, normalized, starting_quote) if state.continuation_out else ""

    normalized_lines = normalize_continuations(lines)
    normalized_lines = normalize_openmp_continuation_sentinels(normalized_lines)
    line_count_before_returns = len(normalized_lines)
    normalized_lines = remove_terminal_procedure_returns(normalized_lines, procedure_cases)
    returns_removed = len(normalized_lines) != line_count_before_returns

    if wrap:
        # RETURN removal is the only operation since the last scope analysis
        # that can shift physical line numbers. Rebuild only when it happened.
        if returns_removed:
            post_return_source = "".join(normalized_lines)
            procedure_cases = extract_procedure_cases(post_return_source)
            scoped_declared_names = extract_scoped_declared_names(post_return_source)
        formatted_lines = rewrap_lines(
            normalized_lines,
            procedure_cases=procedure_cases,
            preserved_names=symbol_cases,
            preprocessor_names=preprocessor_names,
            scoped_declared_names=scoped_declared_names,
            uppercase_single_l=uppercase_single_l,
        )
    else:
        formatted_lines = normalized_lines
    formatted_lines = normalize_declaration_separator_alignment(formatted_lines)
    formatted_lines = normalize_program_unit_spacing(formatted_lines)
    formatted_lines = limit_blank_lines(formatted_lines)
    normalized_output = normalize_output_whitespace("".join(formatted_lines))
    return normalized_output.replace("\n", line_ending)


def lowercase_file(
    path: Path,
    wrap: bool = True,
    module_cases: Mapping[str, str] | None = None,
    symbol_cases: Mapping[str, str] | None = None,
    procedure_cases: Iterable[ProcedureDeclarationCases] = (),
    scope_cases: Iterable[NamedScopeCase] = (),
    type_procedure_cases: Mapping[str, str] | None = None,
    type_component_cases: Mapping[tuple[str, str], str] | None = None,
    variable_type_cases: Mapping[str, str] | None = None,
    type_component_type_cases: Mapping[tuple[str, str], str] | None = None,
    uppercase_single_l: bool = False,
    macro_cases: Mapping[str, str] | Collection[str] = (),
) -> bool:
    """Rewrite *path* and return whether its contents changed."""
    original = read_text_preserving_newlines(path)
    updated = format_text(
        original,
        wrap=wrap,
        module_cases=module_cases,
        symbol_cases=symbol_cases,
        procedure_cases=procedure_cases,
        scope_cases=scope_cases,
        type_procedure_cases=type_procedure_cases,
        type_component_cases=type_component_cases,
        variable_type_cases=variable_type_cases,
        type_component_type_cases=type_component_type_cases,
        uppercase_single_l=uppercase_single_l,
        macro_cases=macro_cases,
    )
    if updated == original:
        return False
    write_text_atomic(path, updated)
    return True


def write_text_atomic(path: Path, text: str) -> None:
    """Atomically replace file contents while preserving a symlink path itself."""
    write_path = path.resolve(strict=True) if path.is_symlink() else path
    mode = write_path.stat().st_mode & 0o7777
    temporary_path: str | None = None
    try:
        file_descriptor, temporary_path = tempfile.mkstemp(prefix=f".{write_path.name}.", dir=write_path.parent)
        with os.fdopen(file_descriptor, "w", encoding="utf-8", newline="") as temporary_file:
            os.chmod(temporary_path, mode)
            temporary_file.write(text)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, write_path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            os.unlink(temporary_path)


def read_text_preserving_newlines(path: Path) -> str:
    """Read a text file without translating its line endings."""
    with path.open("r", encoding="utf-8", newline="") as source_file:
        return source_file.read()


def tracked_fortran_paths() -> list[Path]:
    """Return all tracked free-form Fortran files, including submodules."""
    if REPOSITORY_ROOT is None:
        raise RuntimeError("tracked Fortran paths require a valid Git checkout")
    result = subprocess.run(
        ["git", "ls-files", "--recurse-submodules", "-z", "--", *FORTRAN_SOURCE_PATTERNS],
        cwd=REPOSITORY_ROOT,
        env=_git_env(),
        check=True,
        capture_output=True,
    )
    paths = (Path(path) for path in result.stdout.decode().split("\0") if path)
    return sorted(REPOSITORY_ROOT / path for path in paths)


def read_source_files(paths: Iterable[Path]) -> dict[Path, str]:
    """Read all source files needed for declaration-case resolution."""
    return {path: read_text_preserving_newlines(path) for path in paths}


def macro_name_argument(value: str) -> str:
    """Argparse converter for -D/--define NAME[=VALUE]."""
    name = value.split("=", 1)[0]
    if re.fullmatch(IDENTIFIER, name) is None:
        raise argparse.ArgumentTypeError(f"invalid macro name: {name!r}")
    return name


def _validated_fortran_path(path: Path) -> Path:
    if not path.is_absolute() and REPOSITORY_ROOT is not None and not path.exists():
        candidate = REPOSITORY_ROOT / path
        if candidate.exists():
            path = candidate
    resolved = path.resolve()
    if resolved.suffix.lower() not in FORTRAN_SOURCE_EXTENSIONS:
        extensions = ", ".join(sorted(FORTRAN_SOURCE_EXTENSIONS))
        raise ValueError(f"Expected a free-form Fortran source ({extensions}): {resolved}")
    if not resolved.is_file():
        raise ValueError(f"Fortran source file does not exist: {resolved}")
    return resolved


def _display_path(path: Path) -> Path:
    if REPOSITORY_ROOT is not None:
        try:
            return path.relative_to(REPOSITORY_ROOT)
        except ValueError:
            pass
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path, help="free-form Fortran source files to update in place")
    parser.add_argument(
        "--all", action="store_true", help="update every tracked free-form Fortran source in this repository"
    )
    parser.add_argument(
        "--stdin", action="store_true", help="read Fortran source from standard input and write it to standard output"
    )
    parser.add_argument(
        "--stdout", action="store_true", help="write one input file's formatted source to standard output"
    )
    parser.add_argument(
        "--isolated",
        action="store_true",
        help="resolve declaration cases from the specified files only, without scanning repository sources",
    )
    parser.add_argument("--check", action="store_true", help="report files that would change without rewriting them")
    parser.add_argument("--diff", action="store_true", help="show unified diffs without rewriting files")
    parser.add_argument(
        "--wrap",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="re-wrap code lines longer than 120 columns (default: enabled)",
    )
    parser.add_argument(
        "--uppercase-single-l",
        action="store_true",
        help="replace standalone l/L identifiers with L",
    )
    parser.add_argument(
        "-D",
        "--define",
        dest="macro_names",
        action="append",
        default=[],
        type=macro_name_argument,
        metavar="NAME[=VALUE]",
        help="preserve/use this exact case for a macro supplied externally; may be repeated",
    )
    args = parser.parse_args()
    if args.stdin and (args.paths or args.all or args.stdout or args.check or args.diff or args.isolated):
        parser.error("--stdin cannot be combined with paths, --all, --stdout, --check, --diff, or --isolated")
    if args.stdout and (len(args.paths) != 1 or args.all or args.check or args.diff):
        parser.error("--stdout requires exactly one path and cannot be combined with --all, --check, or --diff")
    if args.all and args.paths:
        parser.error("--all cannot be combined with paths")
    if args.isolated and (args.all or not args.paths):
        parser.error("--isolated requires one or more explicit paths")
    return args


def main() -> int:
    args = parse_args()
    macro_cases = normalize_macro_cases(args.macro_names)
    if args.stdin:
        source = sys.stdin.read()
        cases = collect_declaration_cases({Path("<stdin>"): source})[Path("<stdin>")]
        sys.stdout.write(
            format_text(
                source,
                wrap=args.wrap,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
                type_procedure_cases=cases.type_procedure_cases,
                type_component_cases=cases.type_component_cases,
                variable_type_cases=cases.variable_type_cases,
                type_component_type_cases=cases.type_component_type_cases,
                uppercase_single_l=args.uppercase_single_l,
                macro_cases=macro_cases,
            )
        )
        return 0

    if args.stdout:
        try:
            path = _validated_fortran_path(args.paths[0])
        except ValueError as error:
            print(error, file=sys.stderr)
            return 2
        source_paths = [path]
        if REPOSITORY_ROOT is not None and not args.isolated:
            source_paths = [*tracked_fortran_paths(), path]
        sources = read_source_files(dict.fromkeys(source_paths))
        source = sources[path]
        cases = collect_declaration_cases(sources, target_paths=(path,))[path]
        sys.stdout.write(
            format_text(
                source,
                wrap=args.wrap,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
                type_procedure_cases=cases.type_procedure_cases,
                type_component_cases=cases.type_component_cases,
                variable_type_cases=cases.variable_type_cases,
                type_component_type_cases=cases.type_component_type_cases,
                uppercase_single_l=args.uppercase_single_l,
                macro_cases=macro_cases,
            )
        )
        return 0

    if args.all or not args.paths:
        if REPOSITORY_ROOT is None:
            print("--all or no paths requires a valid Git checkout", file=sys.stderr)
            return 2
        target_paths = tracked_fortran_paths()
    else:
        try:
            target_paths = [_validated_fortran_path(path) for path in args.paths]
        except ValueError as error:
            print(error, file=sys.stderr)
            return 2
    changed_paths: list[Path] = []
    source_paths = (
        [*tracked_fortran_paths(), *target_paths]
        if REPOSITORY_ROOT is not None and not args.isolated and args.paths
        else target_paths
    )
    originals = read_source_files(dict.fromkeys(source_paths))
    declaration_cases = collect_declaration_cases(originals, target_paths=target_paths)

    for path in target_paths:
        original = originals[path]
        cases = declaration_cases[path]
        formatted = format_text(
            original,
            wrap=args.wrap,
            module_cases=cases.module_cases,
            symbol_cases=cases.symbol_cases,
            procedure_cases=cases.procedure_cases,
            scope_cases=cases.scope_cases,
            type_procedure_cases=cases.type_procedure_cases,
            type_component_cases=cases.type_component_cases,
            variable_type_cases=cases.variable_type_cases,
            type_component_type_cases=cases.type_component_type_cases,
            uppercase_single_l=args.uppercase_single_l,
            macro_cases=macro_cases,
        )
        if args.check or args.diff:
            if formatted != original:
                changed_paths.append(path)
                if args.diff:
                    display_path = _display_path(path)
                    sys.stdout.writelines(
                        difflib.unified_diff(
                            original.splitlines(keepends=True),
                            formatted.splitlines(keepends=True),
                            fromfile=f"a/{display_path}",
                            tofile=f"b/{display_path}",
                        )
                    )
        elif formatted != original:
            write_text_atomic(path, formatted)
            changed_paths.append(path)

    for path in changed_paths:
        if not args.diff:
            print(path)
    return int((args.check or args.diff) and bool(changed_paths))


if __name__ == "__main__":
    raise SystemExit(main())
