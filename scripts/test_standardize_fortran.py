import unittest
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import scripts.standardize_fortran as formatter
from scripts.standardize_fortran import (
    _declared_variable_names,
    _iter_code_statements_with_lines,
    collect_declaration_cases,
    extract_declared_names,
    extract_module_variable_names,
    extract_module_variable_types,
    extract_procedure_cases,
    format_text,
    join_lexical_token_continuations,
    lowercase_file,
    parse_args,
    read_source_files,
)


class CommandLineTests(unittest.TestCase):
    def test_invalid_flag_combinations_use_argparse_errors(self) -> None:
        invalid_arguments = (
            (["--stdin", "input.f90"], "--stdin cannot be combined"),
            (["--stdout"], "--stdout requires exactly one path"),
            (["--stdout", "input.f90", "--diff"], "--stdout requires exactly one path"),
            (["--all", "input.f90"], "--all cannot be combined with paths"),
        )
        for arguments, message in invalid_arguments:
            with self.subTest(arguments=arguments):
                stderr = StringIO()
                with (
                    patch("sys.argv", ["standardize_fortran.py", *arguments]),
                    redirect_stderr(stderr),
                    self.assertRaises(SystemExit) as error,
                ):
                    parse_args()
                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())
                self.assertNotIn("Traceback", stderr.getvalue())

    def test_uppercase_single_l_option(self) -> None:
        with patch("sys.argv", ["standardize_fortran.py", "--uppercase-single-l", "input.f90"]):
            self.assertTrue(parse_args().uppercase_single_l)

    def test_isolated_option(self) -> None:
        with patch("sys.argv", ["standardize_fortran.py", "--isolated", "input.f90"]):
            self.assertTrue(parse_args().isolated)

    def test_explicit_path_does_not_require_git_checkout(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.f90"
            path.write_text("use mixed\n")
            with (
                patch.object(formatter, "REPOSITORY_ROOT", None),
                patch.object(formatter, "tracked_fortran_paths", side_effect=AssertionError("unexpected Git lookup")),
                patch("sys.argv", ["standardize_fortran.py", "--check", str(path)]),
            ):
                self.assertEqual(formatter.main(), 0)

    def test_isolated_path_does_not_scan_repository_sources(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.f90"
            path.write_text("x = 1\n")
            with (
                patch.object(formatter, "tracked_fortran_paths", side_effect=AssertionError("unexpected Git lookup")),
                patch("sys.argv", ["standardize_fortran.py", "--isolated", "--check", str(path)]),
            ):
                self.assertEqual(formatter.main(), 0)


class FormattingTests(unittest.TestCase):
    def test_preserves_spacing_in_named_common_blocks(self) -> None:
        source = "      COMMON / RASET1 / U, V\n"
        self.assertEqual(format_text(source, wrap=False), "      common /RASET1/ U, V\n")

    def test_removes_only_redundant_nested_parentheses(self) -> None:
        source = """\
if ((x)) then
    value = (((a + b))) / c
    value = a / ((b + c))
    call work((x))
    value = ((a + b) * c)
    print *, "((literal))"
! ((comment))
#define P ((preprocessor))
end if
"""
        expected = """\
if (x) then
    value = (a + b)/c
    value = a/(b + c)
    call work((x))
    value = ((a + b)*c)
    print *, "((literal))"
! ((comment))
#define P ((preprocessor))
end if
"""
        self.assertEqual(format_text(source, wrap=False), expected)

    def test_normalizes_dimension_and_write_output_spacing(self) -> None:
        source = """\
integer, dimension (:) :: values
write(*, *)'Warning...'
write(unit, '(1I6,4E15.6)')il, value
write(unit, '(1I6,4E15.6)')
write(unit, '(1I6,4E15.6)') &
write(unit, '(1I6,4E15.6)' ) ! no output item
print *, "write(*)'literal'"
! write(*)'comment'
"""
        expected = """\
integer, dimension(:) :: values
write(*, *) 'Warning...'
write(unit, '(1I6,4E15.6)') il, value
write(unit, '(1I6,4E15.6)')
write(unit, '(1I6,4E15.6)') &
write(unit, '(1I6,4E15.6)' ) ! no output item
print *, "write(*)'literal'"
! write(*)'comment'
"""
        self.assertEqual(format_text(source, wrap=False), expected)

    def test_preserves_nested_parentheses_in_arguments_and_associations(self) -> None:
        source = "x = f((x))\nassociate(y => ((x)))\n"
        self.assertEqual(formatter.remove_redundant_nested_parentheses(source), source)


class DeclarationCaseTests(unittest.TestCase):
    def test_conditional_sentinel_body_uses_declared_case(self) -> None:
        source = """\
module t
integer :: MyVar
contains
subroutine s()
!$ myvar = 1
myvar = 2
end subroutine s
end module t
"""
        cases = collect_declaration_cases({Path("sentinel.f90"): source})[Path("sentinel.f90")]
        self.assertEqual(
            format_text(
                source,
                wrap=False,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
            ),
            """\
module t
integer :: MyVar

contains

subroutine s
!$ MyVar = 1
MyVar = 2

end subroutine s

end module t
""",
        )

    def test_declaration_array_constructor_is_one_entity(self) -> None:
        self.assertEqual(_declared_variable_names("integer :: arr(2) = [A, B]"), ["arr"])

    def test_extracts_and_matches_declarations(self) -> None:
        declarations = """\
module MiXeD
    type :: MyType
    function MyFunc(value)
    end function MyFunc
end module MiXeD
"""
        uses = """\
use mixed
type(mytype) :: value
call myfunc(value)
character(len=6) :: text = "mytype"
! myfunc and mytype stay in this comment
"""
        self.assertEqual(
            [(item.kind, item.name) for item in extract_declared_names(declarations)],
            [("module", "MiXeD"), ("type", "MyType"), ("procedure", "MyFunc")],
        )
        cases = collect_declaration_cases({Path("decl.f90"): declarations, Path("use.f90"): uses})[Path("use.f90")]
        self.assertEqual(
            format_text(uses, wrap=False, module_cases=cases.module_cases, symbol_cases=cases.symbol_cases),
            """\
use MiXeD
type(MyType) :: value
call MyFunc(value)
character(len=6) :: text = "mytype"
! myfunc and mytype stay in this comment
""",
        )

    def test_duplicate_resolution(self) -> None:
        sources = {
            Path("one.f90"): "module First\n",
            Path("two.f90"): "module FIRST\n",
            Path("same.f90"): "module Same\nmodule Same\n",
            Path("uses.f90"): "use first\n",
            Path("duplicate.f90"): "module Foo\nmodule fOO\nuse foo\n",
            Path("external.f90"): "module ExternalCase\n",
            Path("current.f90"): "module externalcase\nuse EXTERNALCASE\n",
        }
        cases = collect_declaration_cases(sources)
        self.assertNotIn("first", cases[Path("uses.f90")].module_cases)
        self.assertNotIn("foo", cases[Path("duplicate.f90")].module_cases)
        self.assertEqual(cases[Path("same.f90")].module_cases["same"], "Same")
        self.assertEqual(cases[Path("current.f90")].module_cases["externalcase"], "externalcase")

    def test_named_end_cases_follow_start_cases(self) -> None:
        source = """\
module MiXeD
contains
subroutine Work
end subroutine work
end module mixed
"""
        cases = collect_declaration_cases({Path("scope.f90"): source})[Path("scope.f90")]
        self.assertEqual(
            format_text(
                source,
                wrap=False,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
            ),
            """\
module MiXeD

contains

subroutine Work

end subroutine Work

end module MiXeD
""",
        )

    def test_compact_named_ends_follow_start_cases(self) -> None:
        source = """\
module MiXeD
contains
subroutine Work
endsubroutine work
endmodule mixed
"""
        cases = collect_declaration_cases({Path("scope.f90"): source})[Path("scope.f90")]
        self.assertEqual(
            format_text(
                source,
                wrap=False,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
            ),
            """\
module MiXeD

contains

subroutine Work

end subroutine Work

end module MiXeD
""",
        )

    def test_compact_procedure_end_closes_local_case_scope(self) -> None:
        source = """\
subroutine Work(F)
    real :: f, D
    F = d
endsubroutine work
subroutine Other
    real :: value
    VALUE = 1
endsubroutine other
"""
        cases = extract_procedure_cases(source)
        self.assertEqual([(case.start_line, case.end_line) for case in cases], [(0, 3), (4, 7)])

    def test_nested_procedure_uses_innermost_local_case(self) -> None:
        source = """\
subroutine outer(size)
    integer :: size
contains
    subroutine inner(size)
        integer :: SIZE
        size = 1
    end subroutine inner
end subroutine outer
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
subroutine outer(size)
    integer :: size

contains

    subroutine inner(SIZE)
        integer :: SIZE
        SIZE = 1

    end subroutine inner

end subroutine outer
""",
        )

    def test_bare_end_closes_local_case_scope(self) -> None:
        source = """\
subroutine first
    integer :: Size
    Size = 1
end
subroutine second
    integer :: value
    Size = 2
end
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
subroutine first
    integer :: Size
    Size = 1

end
subroutine second
    integer :: value
    size = 2

end
""",
        )

    def test_module_procedure_is_not_a_module_or_module_variable_scope(self) -> None:
        source = """\
submodule (parent) child
contains
module subroutine Foo(SIZE)
    integer :: SIZE
    SIZE = 1

end subroutine Foo
end submodule child
"""
        self.assertEqual([(item.kind, item.name) for item in extract_declared_names(source)], [("procedure", "Foo")])
        self.assertEqual(extract_procedure_cases(source)[0].local_cases, {"size": "SIZE"})
        self.assertEqual(extract_module_variable_names(source), [])
        self.assertEqual(
            format_text(source, wrap=False),
            """\
submodule (parent) child
contains
module subroutine Foo(SIZE)
    integer :: SIZE
    SIZE = 1

end subroutine Foo
end submodule child
""",
        )

    def test_interface_dummies_are_not_module_variables(self) -> None:
        source = """\
module M
    interface
        subroutine ext(ArgCase)
            type(Widget) :: ArgCase
        end subroutine ext
    end interface
contains
    subroutine s
        print *, argcase
    end subroutine s
end module M
"""
        self.assertEqual(extract_module_variable_names(source), [])
        self.assertEqual(extract_module_variable_types(source), {})
        cases = collect_declaration_cases({Path("interface.f90"): source})[Path("interface.f90")]
        self.assertEqual(
            format_text(source, wrap=False, symbol_cases=cases.symbol_cases),
            """\
module M
    interface
        subroutine ext(ArgCase)
            type(Widget) :: ArgCase
        end subroutine ext
    end interface

contains

    subroutine s
        print *, argcase

    end subroutine s

end module M
""",
        )

    def test_lexical_join_refreshes_procedure_line_ranges(self) -> None:
        source = """\
integer :: ab&
&cd
subroutine foo(SIZE)
    integer :: SIZE
    SIZE=1
end subroutine foo
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
integer :: abcd
subroutine foo(SIZE)
    integer :: SIZE
    SIZE = 1

end subroutine foo
""",
        )

    def test_local_variables_override_and_leave_global_pool(self) -> None:
        sources = {
            Path("global.f90"): "function f()\nend function f\n",
            Path("local.f90"): """\
subroutine Work(F)
    real :: f, D
    type(MyType) legacy
    F = d
    LEGACY = F
end subroutine work
""",
            Path("program.f90"): """\
program Demo
    real :: x
    X = x
end program demo
""",
        }
        cases = collect_declaration_cases(sources)
        local = cases[Path("local.f90")]
        self.assertEqual(local.symbol_cases["f"], "f")
        self.assertEqual(local.procedure_cases[0].local_cases, {"f": "f", "d": "D", "legacy": "legacy"})
        self.assertEqual(
            format_text(
                sources[Path("local.f90")],
                wrap=False,
                module_cases=local.module_cases,
                symbol_cases=local.symbol_cases,
                procedure_cases=local.procedure_cases,
                scope_cases=local.scope_cases,
            ),
            """\
subroutine Work(f)
    real :: f, D
    type(MyType) legacy
    f = D
    legacy = f

end subroutine Work
""",
        )
        program = cases[Path("program.f90")]
        self.assertEqual(program.procedure_cases[0].local_cases, {"x": "x"})
        self.assertEqual(
            format_text(
                sources[Path("program.f90")],
                wrap=False,
                module_cases=program.module_cases,
                symbol_cases=program.symbol_cases,
                procedure_cases=program.procedure_cases,
                scope_cases=program.scope_cases,
            ),
            """\
program Demo
    real :: x
    x = x

end program Demo
""",
        )

    def test_local_case_does_not_apply_to_derived_type_components(self) -> None:
        sources = {
            Path("global.f90"): """\
type :: ComponentCase
end type ComponentCase
""",
            Path("components.f90"): """\
subroutine Work
    real :: WINDOW
    WINDOW = RedWin%componentcase%Window_f_a(a, winamp)
end subroutine work
""",
        }
        cases = collect_declaration_cases(sources)[Path("components.f90")]
        self.assertEqual(
            format_text(
                sources[Path("components.f90")].replace("WINDOW =", "window ="),
                wrap=False,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
            ),
            """\
subroutine Work
    real :: WINDOW
    WINDOW = RedWin%ComponentCase%Window_f_a(a, winamp)

end subroutine Work
""",
        )

    def test_module_variables_are_case_matched_without_leaking_local_shadowing(self) -> None:
        sources = {
            Path("config.f90"): """\
module config
    integer :: FeedbackLevel
    type :: State
        real :: transfer_times(:)
        real :: H0
    end type State
contains
end module config
""",
            Path("uses.f90"): """\
module Uses
    use config

contains

    subroutine Work(Feedbacklevel, H0)
        integer, intent(in) :: Feedbacklevel
        real, intent(in)    :: H0
        type(State)         :: obj, P
        real                :: x
        x = obj%Transfer_Times(1)
        x = P%h0
        if (feedbacklevel > 0) then
        end if
    end subroutine work
end module uses
""",
            Path("external.f90"): """\
subroutine External
    if (feedbacklevel > 0) then
    end if
end subroutine external
""",
        }
        cases = collect_declaration_cases(sources)
        self.assertEqual(cases[Path("external.f90")].symbol_cases["feedbacklevel"], "FeedbackLevel")
        self.assertEqual(cases[Path("external.f90")].symbol_cases["transfer_times"], "transfer_times")
        self.assertEqual(
            format_text(
                sources[Path("uses.f90")],
                wrap=False,
                module_cases=cases[Path("uses.f90")].module_cases,
                symbol_cases=cases[Path("uses.f90")].symbol_cases,
                procedure_cases=cases[Path("uses.f90")].procedure_cases,
                scope_cases=cases[Path("uses.f90")].scope_cases,
            ),
            """\
module Uses
    use config

contains

    subroutine Work(Feedbacklevel, H0)
        integer, intent(in) :: Feedbacklevel
        real, intent(in)    :: H0
        type(State)         :: obj, P
        real                :: x
        x = obj%transfer_times(1)
        x = P%H0
        if (Feedbacklevel > 0) then
        end if

    end subroutine Work

end module Uses
""",
        )

    def test_select_type_aliases_are_local(self) -> None:
        sources = {
            Path("global.f90"): "module config\n    integer :: SP\nend module config\n",
            Path("uses.f90"): """\
subroutine Work(obj)
    class(*) :: obj
    select type (sp => obj)
    class is (SomeType)
        call sp%init()
    end select
end subroutine work
""",
        }
        cases = collect_declaration_cases(sources)[Path("uses.f90")]
        self.assertEqual(cases.procedure_cases[0].local_cases["sp"], "sp")
        self.assertEqual(
            format_text(
                sources[Path("uses.f90")],
                wrap=False,
                module_cases=cases.module_cases,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                scope_cases=cases.scope_cases,
            ),
            """\
subroutine Work(obj)
    class(*) :: obj
    select type (sp => obj)
    class is (SomeType)
        call sp%init()
    end select

end subroutine Work
""",
        )


class ContinuationTests(unittest.TestCase):
    def test_file_reads_and_writes_preserve_crlf(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.f90"
            path.write_bytes(b"X=1  \r\n")

            self.assertEqual(read_source_files([path])[path], "X=1  \r\n")
            self.assertTrue(lowercase_file(path, wrap=False))
            self.assertEqual(path.read_bytes(), b"X = 1\r\n")

    def test_preserves_continued_cpp_directives(self) -> None:
        source = """\
#define MAX 42
#define X \\
MAX \\
VALUE
IF (X .EQ. MAX) THEN
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
#define MAX 42
#define X \\
MAX \\
VALUE
if (X == MAX) then
""",
        )

    def test_joins_continuations_inside_lexical_tokens(self) -> None:
        source = """\
program p
    long&
&identifier = 1
end program p
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
program p
    longidentifier = 1

end program p
""",
        )

    def test_continuation_whitespace_separates_tokens(self) -> None:
        source = """\
subrou&
&tine Work
call &
&work()
"""
        self.assertEqual(
            join_lexical_token_continuations(source),
            """\
subroutine Work
call &
&work()
""",
        )
        self.assertEqual(
            [statement.text for statement in _iter_code_statements_with_lines(source)],
            ["subroutine Work", "call work()"],
        )

    def test_preserves_lexical_continuation_markers_around_inline_comment(self) -> None:
        source = """\
subrou& ! split token
&tine s()

end subroutine s
"""
        self.assertEqual(join_lexical_token_continuations(source), source)
        self.assertEqual(format_text(source, wrap=False), source)

    def test_leaves_continued_character_literals_unchanged(self) -> None:
        source = 'value = "a&\n&b"\n'
        self.assertEqual(format_text(source, wrap=False), source)

    def test_splits_top_level_semicolon_statements(self) -> None:
        source = "x = 1; call work('a;b'); y = (2; 3)\n"
        self.assertEqual(
            [statement.text for statement in _iter_code_statements_with_lines(source)],
            ["x = 1", "call work(#####)", "y = (2; 3)"],
        )


class SpacingTests(unittest.TestCase):
    def test_intrinsic_names_do_not_override_local_variables(self) -> None:
        source = """\
subroutine Work
    integer :: Size
    Size = 1
    print *, Size
    value = SIN(x) + MAX(x, 1)
end subroutine Work
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
subroutine Work
    integer :: Size
    Size = 1
    print *, Size
    value = sin(x) + max(x, 1)

end subroutine Work
""",
        )

    def test_global_symbols_do_not_override_intrinsics_or_real_exponents(self) -> None:
        source = """\
n = SIZE(values)
real(dl), parameter :: eps_xf = 5E-4 + 1.E-5_dl + 1.D-6
"""
        self.assertEqual(
            format_text(source, wrap=False, symbol_cases={"size": "Size", "e": "E"}),
            """\
n = size(values)
real(dl), parameter :: eps_xf = 5e-4 + 1.e-5_dl + 1.d-6
""",
        )

    def test_module_declared_names_keep_their_case_against_specifiers_and_intrinsics(self) -> None:
        source = """\
module M
    integer :: SIZE
    integer :: MAX
contains
    subroutine S
        print *, SIZE
        value = MAX(a, b)
        inquire (unit=1, size=SIZE)
    end subroutine S
end module M
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
module M
    integer :: SIZE
    integer :: MAX

contains

    subroutine S
        print *, SIZE
        value = MAX(a, b)
        inquire(unit=1, size=SIZE)

    end subroutine S

end module M
""",
        )

    def test_declared_names_do_not_leak_across_modules_in_the_same_file(self) -> None:
        source = """\
module A
    integer, private :: SIZE
contains
    subroutine set_size(n)
        integer, intent(in) :: n
        SIZE = n
    end subroutine set_size
end module A

module B
contains
    subroutine report(x)
        real, intent(in) :: x(:)
        print *, SIZE(x)
    end subroutine report
end module B
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
module A
    integer, private :: SIZE

contains

    subroutine set_size(n)
        integer, intent(in) :: n
        SIZE = n

    end subroutine set_size

end module A

module B

contains

    subroutine report(x)
        real, intent(in) :: x(:)
        print *, size(x)

    end subroutine report

end module B
""",
        )

    def test_declared_names_do_not_leak_from_a_type_component(self) -> None:
        source = """\
module C
    type Foo
        integer :: SIZE
    end type Foo
contains
    subroutine report(x)
        real, intent(in) :: x(:)
        print *, SIZE(x)
    end subroutine report
end module C
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
module C
    type Foo
        integer :: SIZE
    end type Foo

contains

    subroutine report(x)
        real, intent(in) :: x(:)
        print *, size(x)

    end subroutine report

end module C
""",
        )

    def test_component_cases_apply_inside_modules(self) -> None:
        sources = {
            Path("container.f90"): """\
module demo
    type :: container
        integer :: CAMB_PK
    end type container
end module demo
""",
            Path("use_container.f90"): """\
module use_demo

contains
    subroutine clear(value)
        type(container) :: value
        value%CAMB_Pk = 0
    end subroutine clear
end module use_demo
""",
        }
        cases = collect_declaration_cases(sources)[Path("use_container.f90")]
        self.assertEqual(
            format_text(
                sources[Path("use_container.f90")],
                wrap=False,
                type_component_cases=cases.type_component_cases,
            ),
            """\
module use_demo

contains

    subroutine clear(value)
        type(container) :: value
        value%CAMB_PK = 0

    end subroutine clear

end module use_demo
""",
        )

    def test_lowercases_standard_statement_specifiers(self) -> None:
        source = """\
module demo


contains
subroutine work(values, prototype, status, message)
    integer :: values, prototype, status
    character(*) :: message
    allocate(values, SOURCE=prototype, MOLD=prototype, STAT=status, ERRMSG=message)
    inquire (UNIT=1, ACCESS='stream', NAME=message)
    read (UNIT=1, REC=1, IOSTAT=status)
    stop 1, QUIET=.TRUE.
end subroutine work
end module demo
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
module demo

contains

subroutine work(values, prototype, status, message)
    integer :: values, prototype, status
    character(*) :: message
    allocate(values, source=prototype, mold=prototype, stat=status, errmsg=message)
    inquire(unit=1, access='stream', name=message)
    read(unit=1, rec=1, iostat=status)
    stop 1, quiet = .true.

end subroutine work

end module demo
""",
        )

    def test_preserves_concatenation_spacing_on_a_continuation_line(self) -> None:
        source = """\
call MpiStop('SP(k) cannot be combined with HMCode_A_baryon/' &
    // 'HMCode_eta_baryon baryonic corrections in HMCode 2015/2016')
"""
        self.assertEqual(format_text(source, wrap=False), source)

    def test_parenthesized_statements_lowercase_unless_locally_shadowed(self) -> None:
        source = """\
WRITE (*, *) value
READ (unit, *) value
OPEN (newunit=unit, file=name)
BACKSPACE (unit)
ALLOCATED (value)
C%Write (*, *) value
subroutine s
    procedure :: Write
    call WRITE()
end subroutine s
"""
        self.assertEqual(
            format_text(source, wrap=False, symbol_cases={"write": "Write"}),
            """\
write(*, *) value
read(unit, *) value
open(newunit=unit, file=name)
backspace(unit)
allocated(value)
C%Write(*, *) value
subroutine s
    procedure :: Write
    call Write()

end subroutine s
""",
        )

    def test_type_bound_procedures_only_supply_component_case(self) -> None:
        sources = {
            Path("type.f90"): """\
module type_module
    type :: State
    contains
        procedure :: BuildValue
    end type State
end module type_module
""",
            Path("use.f90"): """\
call buildvalue()
state%buildvalue()
""",
        }
        cases = collect_declaration_cases(sources)[Path("use.f90")]
        self.assertNotIn("buildvalue", cases.symbol_cases)
        self.assertEqual(cases.type_procedure_cases, {("state", "buildvalue"): "BuildValue"})
        self.assertEqual(
            format_text(
                sources[Path("use.f90")],
                wrap=False,
                symbol_cases=cases.symbol_cases,
                type_procedure_cases=cases.type_procedure_cases,
            ),
            """\
call buildvalue()
State%buildvalue()
""",
        )

    def test_type_bound_procedure_case_uses_the_governing_owner(self) -> None:
        sources = {
            Path("types.f90"): """\
module types
    type :: ThermoData
    contains
        procedure :: values
    end type ThermoData
    type :: OtherData
    contains
        procedure :: Values
    end type OtherData
end module types
""",
            Path("use.f90"): """\
type(ThermoData) :: data
call data%Values()
""",
        }
        cases = collect_declaration_cases(sources)[Path("use.f90")]
        self.assertEqual(
            format_text(
                sources[Path("use.f90")],
                wrap=False,
                type_procedure_cases=cases.type_procedure_cases,
                variable_type_cases=cases.variable_type_cases,
                type_component_type_cases=cases.type_component_type_cases,
            ),
            """\
type(ThermoData) :: data
call data%values()
""",
        )

    def test_old_style_typed_local_entities_govern_case(self) -> None:
        self.assertEqual(_declared_variable_names("type(EvolutionVars) :: EVOut"), ["EVOut"])
        sources = {
            Path("other.f90"): "module other\n    real :: PK\nend module other\n",
            Path("use.f90"): """\
subroutine load
    real(dl) kh, Pk
    read *, PK
end subroutine load
""",
        }
        cases = collect_declaration_cases(sources)[Path("use.f90")]
        self.assertEqual(
            format_text(
                sources[Path("use.f90")],
                wrap=False,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
            ),
            """\
subroutine load
    real(dl) kh, Pk
    read *, Pk

end subroutine load
""",
        )

    def test_top_level_parameter_governs_file_case(self) -> None:
        sources = {
            Path("other.f90"): "module other\n    integer, parameter :: BJL_recurrence_MAX_L = 25\nend module other\n",
            Path("use.f90"): """\
integer, parameter :: BJL_RECURRENCE_MAX_L = 25
if (l > bjl_recurrence_max_l) then
end if
""",
        }
        cases = collect_declaration_cases(sources)[Path("use.f90")]
        self.assertEqual(
            format_text(
                sources[Path("use.f90")],
                wrap=False,
                symbol_cases=cases.symbol_cases,
            ),
            """\
integer, parameter :: BJL_RECURRENCE_MAX_L = 25
if (l > BJL_RECURRENCE_MAX_L) then
end if
""",
        )

    def test_declaration_entities_are_not_replaced_by_global_symbol_case(self) -> None:
        self.assertEqual(
            format_text(
                "type :: T\ncontains\n    procedure :: Error\nend type T\n", wrap=False, symbol_cases={"error": "error"}
            ),
            "type :: T\ncontains\n    procedure :: Error\nend type T\n",
        )

    def test_normalizes_control_keywords_and_bracket_spacing(self) -> None:
        source = """\
endif
elseif  ( x )
if  ( x )then
end   subroutine Work
do   i=1, 3
do  while  ( x )
associate ( x => y )
type (State) :: state
class (Base) :: base
real function value() result ( answer )
if (any(x == y))call work()
read (unit, *) value
write (*, *) value
open (newunit=unit, file=name)
close (unit)
allocate (values(n))
deallocate (values)
inquire (unit=unit, opened=opened)
if (allocated (values)) then
use thing, only : x, y
iqg = - Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
call work( x, [ y ] )
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
end if
else if (x)
if (x) then

end subroutine Work
do i = 1, 3
do while (x)
associate(x => y)
type(State) :: state
class(Base) :: base
real function value() result(answer)
if (any(x == y)) call work()
read(unit, *) value
write(*, *) value
open(newunit=unit, file=name)
close(unit)
allocate(values(n))
deallocate(values)
inquire(unit=unit, opened=opened)
if (allocated(values)) then
use thing, only: x, y
iqg = -Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
call work(x, [y])
""",
        )
        self.assertEqual(format_text("end   function   state_function\n", wrap=False), "end function state_function\n")

    def test_removes_empty_subroutine_arguments_and_spaces_select_type(self) -> None:
        source = """\
subroutine first()
    select type(value)
    end select
end subroutine first
real function second() result(answer)
    select type is (value)
    end select
end function second
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
subroutine first
    select type (value)
    end select

end subroutine first
real function second() result(answer)
    select type is (value)
    end select

end function second
""",
        )

    def test_limits_blank_lines_inside_module_interfaces(self) -> None:
        source = """\
module demo


interface


end interface


contains
subroutine work
end subroutine work
end module demo
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
module demo

interface

end interface

contains

subroutine work

end subroutine work

end module demo
""",
        )

    def test_keeps_exactly_one_blank_line_around_contains(self) -> None:
        source = """\
module demo
integer :: value



contains



subroutine work
end subroutine work
end module demo
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
module demo
integer :: value

contains

subroutine work

end subroutine work

end module demo
""",
        )

    def test_keeps_blank_line_after_contains_following_select_type(self) -> None:
        source = """\
function format_value(value) result(text)
    class(*) :: value
    character(:), allocatable :: text
    select type (value)
    type is (integer)
        text = 'integer'
    end select
    contains
    subroutine error
    end subroutine error
end function format_value
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
function format_value(value) result(text)
    class(*) :: value
    character(:), allocatable :: text
    select type (value)
    type is (integer)
        text = 'integer'
    end select

    contains

    subroutine error

    end subroutine error

end function format_value
""",
        )

    def test_resolves_chained_component_cases(self) -> None:
        sources = {
            Path("types.f90"): """\
module types
    type :: TRanges
        integer :: Count
    end type TRanges
    type :: CAMBdata
        type(TRanges) :: TimeSteps
    end type CAMBdata
    class(CAMBdata), pointer :: State
end module types
""",
            Path("use_types.f90"): """\
module use_types
    use types
contains
    subroutine work
        State%TimeSteps%count = 1
    end subroutine work
end module use_types
""",
        }
        cases = collect_declaration_cases(sources)[Path("use_types.f90")]
        self.assertEqual(
            format_text(
                sources[Path("use_types.f90")],
                wrap=False,
                symbol_cases=cases.symbol_cases,
                type_component_cases=cases.type_component_cases,
                variable_type_cases=cases.variable_type_cases,
                type_component_type_cases=cases.type_component_type_cases,
            ),
            """\
module use_types
    use types

contains

    subroutine work
        State%TimeSteps%Count = 1

    end subroutine work

end module use_types
""",
        )

    def test_resolves_component_case_after_a_local_object_chain(self) -> None:
        sources = {
            Path("types.f90"): """\
module types
    type :: Params
        real :: TCMB
    end type Params
    type :: Data
        type(Params) :: CP
    end type Data
end module types
""",
            Path("use_types.f90"): """\
module use_types
    use types
contains
    subroutine work(this)
        type(Data) :: this
        value = this%CP%tcmb
    end subroutine work
end module use_types
""",
        }
        cases = collect_declaration_cases(sources)[Path("use_types.f90")]
        self.assertEqual(
            format_text(
                sources[Path("use_types.f90")],
                wrap=False,
                symbol_cases=cases.symbol_cases,
                procedure_cases=cases.procedure_cases,
                type_component_cases=cases.type_component_cases,
                variable_type_cases=cases.variable_type_cases,
                type_component_type_cases=cases.type_component_type_cases,
            ),
            """\
module use_types
    use types

contains

    subroutine work(this)
        type(data) :: this
        value = this%CP%TCMB

    end subroutine work

end module use_types
""",
        )

    def test_normalizes_trailing_whitespace_and_file_endings(self) -> None:
        for newline in ("\n", "\r\n", "\r"):
            with self.subTest(newline=repr(newline)):
                source = f"x = 1  {newline}{newline}{newline}"
                self.assertEqual(format_text(source, wrap=False), f"x = 1{newline}")
        self.assertEqual(format_text("x = 1  ", wrap=False), "x = 1\n")
        self.assertEqual(format_text("x = 1\n\n", wrap=False), "x = 1\n")

    def test_does_not_change_spacing_inside_literals_or_comments(self) -> None:
        source = 'text = "endif  ( x ) associate ( y ) only :" ! endif  ( x )\n'
        self.assertEqual(format_text(source, wrap=False), source)

    def test_only_formats_comments_that_start_with_assignment(self) -> None:
        source = """\
! C++ example: x=1
! An explanation with x = 1
! value=first+second
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
! C++ example: x=1
! An explanation with x = 1
! value = first + second
""",
        )

    def test_compound_keywords_are_only_expanded_at_statement_start(self) -> None:
        source = """\
value = endif + enddo + endmodule + elseif
call work(endif, enddo, endmodule, elseif)
endif
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
value = endif + enddo + endmodule + elseif
call work(endif, enddo, endmodule, elseif)
end if
""",
        )

    def test_normalizes_go_to_to_goto(self) -> None:
        source = """\
GO TO 10
if (x) Go To 20
print *, "GO TO 30" ! GO TO 40
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
goto 10
if (x) goto 20
print *, "GO TO 30" ! GO TO 40
""",
        )

    def test_normalizes_post_f2008_language_keywords(self) -> None:
        source = """\
IMPURE  ELEMENTAL FUNCTION f(x)
PURE   ELEMENTAL SUBROUTINE s
CONTIGUOUS :: x
CRITICAL(STAT = istat)
CHANGE   TEAM(newteam)
SELECT  RANK(a)
RANK  DEFAULT
FORM  TEAM(n, team, STAT=istat)
SYNC  ALL(STAT=istat)
SYNC   TEAM(team)
EVENT  POST(event)
EVENT WAIT(event, UNTIL_COUNT =n)
FAIL  IMAGE
LOCK(lockvar, ACQUIRED_LOCK = acquired)
UNLOCK(lockvar)
DO  CONCURRENT(i=1:n) LOCAL_INIT(x) SHARED(y) REDUCE(+:z)
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
impure elemental function f(x)
pure elemental subroutine s
contiguous :: x
critical(stat=istat)
change team (newteam)
select rank (a)
rank default
form team (n, team, stat=istat)
sync all(stat=istat)
sync team (team)
event post(event)
event wait(event, until_count=n)
fail image
lock(lockvar, acquired_lock=acquired)
unlock(lockvar)
do concurrent(i=1:n) local_init(x) shared(y) reduce(+:z)
""",
        )

    def test_normalizes_not_operator_and_comment_spacing(self) -> None:
        source = """\
if(.not.x.and.y.or.z)then!inline comment
!full comment
!   already a comment
call work()
!   single comment
call work()
!$omp parallel do,private(i,j)
!$omp do default(shared) schedule(dynamic)
!$ use omp_lib,only: omp_get_thread_num
!$acc parallel loop private(i,j)
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
if (.not. x .and. y .or. z) then ! inline comment
! full comment
!   already a comment
call work()
! single comment
call work()
!$OMP PARALLEL DO, PRIVATE(i, j)
!$OMP DO DEFAULT(SHARED), SCHEDULE(DYNAMIC)
!$ use omp_lib, only: omp_get_thread_num
!$acc parallel loop private(i, j)
""",
        )

    def test_lowercases_logical_operators(self) -> None:
        source = "if (.NOT. x .AND. y .OR. z) call work()\n"
        self.assertEqual(format_text(source, wrap=False), "if (.not. x .and. y .or. z) call work()\n")

    def test_preserves_case_sensitive_preprocessor_macros(self) -> None:
        source = """\
#define SIZE 4
program p
  implicit none
  integer :: a(SIZE)
  a = 1
  print *, SIZE

end program p
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
#define SIZE 4
program p
  implicit none
  integer :: a(SIZE)
  a = 1
  print *, SIZE

end program p
""",
        )
        self.assertEqual(
            format_text("#define VALUE 4\nprint *, VALUE\n", wrap=False, symbol_cases={"value": "Value"}),
            "#define VALUE 4\nprint *, VALUE\n",
        )

    def test_optionally_uppercases_single_l_and_modernizes_array_constructors(self) -> None:
        source = "l = (/ 1, 2, 3 /)\nprint *, l\n"
        self.assertEqual(format_text(source, wrap=False), "l = [1, 2, 3]\nprint *, l\n")
        self.assertEqual(format_text(source, wrap=False, uppercase_single_l=True), "L = [1, 2, 3]\nprint *, L\n")

    def test_modernizes_multiline_array_constructors(self) -> None:
        source = """\
real :: values(2) = (/ &
  1.0, &
  2.0 /)
10 format(2x, &
 /)
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
real :: values(2) = [ &
  1.0, &
  2.0]
10 format(2x, &
 /)
""",
        )

    def test_removes_terminal_function_and_subroutine_returns(self) -> None:
        source = """\
subroutine work
    return
end
function value
    value = 1
    RETURN
end function value
program example
    return
end program example
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
subroutine work

end
function value
    value = 1

end function value
program example
    return

end program example
""",
        )

    def test_spaces_bare_program_unit_ends_like_named_ends(self) -> None:
        source = """\
subroutine first
    integer :: value
    value = 1
end
subroutine second
    integer :: value
    value = 2
end
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
subroutine first
    integer :: value
    value = 1

end
subroutine second
    integer :: value
    value = 2

end
""",
        )

    def test_preserves_compiler_directive_sentinels_without_comment_spacing(self) -> None:
        source = """\
!$acc parallel
!DIR$ IVDEP
!DEC$ ATTRIBUTES FORCEINLINE :: Work
!GCC$ ATTRIBUTES HOT :: Work
! ordinary comment
"""
        self.assertEqual(format_text(source, wrap=False), source)

    def test_wraps_openmp_directives(self) -> None:
        source = """\
                    !$omp parallel do default(shared) schedule(static) private(Cl, ell, reall, fac, n, LimbRec, LimbRec2)
"""
        formatted = format_text(source)
        self.assertEqual(
            formatted,
            """\
                    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), &
                    !$OMP PRIVATE(Cl, ell, reall, fac, n, LimbRec, LimbRec2)
""",
        )
        self.assertEqual(format_text(formatted), formatted)

    def test_preserves_openmp_breaks_and_canonicalizes_continuation_sentinels(self) -> None:
        source = """\
!$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC), IF(n >= MATHUTILS_OMP_GAUSS_THRESHOLD) &
!$OMP& PRIVATE(i, j, k, iter, p1, p2, p3, pp, z, dz, wi) &
!$OMP & SHARED(x, w, n, m, failed)
"""
        self.assertEqual(
            format_text(source),
            """\
!$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC), IF(n >= MATHUTILS_OMP_GAUSS_THRESHOLD) &
!$OMP PRIVATE(i, j, k, iter, p1, p2, p3, pp, z, dz, wi) &
!$OMP SHARED(x, w, n, m, failed)
""",
        )

    def test_limits_blank_lines_and_preserves_declaration_alignment(self) -> None:
        source = """\
real      :: first
integer   :: second

! A comment does not interrupt declaration alignment.

logical   :: third



real   :: unaligned

real ::   padded_one
integer ::   padded_two
character ::    unaligned_padding
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
real    :: first
integer :: second

! A comment does not interrupt declaration alignment.

logical :: third


real :: unaligned

real :: padded_one
integer :: padded_two
character :: unaligned_padding
""",
        )

    def test_normalizes_old_style_declaration_spacing_and_optional_order(self) -> None:
        source = """\
    real(dl)  x
    real(dl)kh, PK
    real(dp), optional, intent(out) :: sin_k
    real(dp), intent(in), optional :: cos_k
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
    real(dl) x
    real(dl) kh, PK
    real(dp), intent(out), optional :: sin_k
    real(dp), intent(in), optional :: cos_k
""",
        )

    def test_reduces_declaration_alignment_to_its_minimum(self) -> None:
        source = """\
    procedure, private  :: WriteSizedArray1
    procedure, private  :: WriteSizedArray2
    procedure, private  :: ReadSizedArray_R

    generic  :: LoadTxt => LoadTxt_2D, LoadTxt_1D
    generic  :: SaveTxt => WriteTextMatrix, WriteTextVector

    integer, intent(in)   :: md
    integer, intent(in)   :: nxd
    real(GI), intent(in)    :: xd(nxd)
    real(GI), intent(in)  :: zd(nxd)

    real(GI), intent(in)  :: xx1
    real(GI)              :: fn_val
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
    procedure, private :: WriteSizedArray1
    procedure, private :: WriteSizedArray2
    procedure, private :: ReadSizedArray_R

    generic :: LoadTxt => LoadTxt_2D, LoadTxt_1D
    generic :: SaveTxt => WriteTextMatrix, WriteTextVector

    integer, intent(in)  :: md
    integer, intent(in)  :: nxd
    real(GI), intent(in) :: xd(nxd)
    real(GI), intent(in) :: zd(nxd)

    real(GI), intent(in) :: xx1
    real(GI)             :: fn_val
""",
        )

    def test_alignment_compresses_through_comment_lines(self) -> None:
        source = """\
    real(dl), intent(in)              :: ax      !! left endpoint of initial interval
    real(dl), intent(in)              :: bx      !! right endpoint of initial interval
    real(dl), intent(in)              :: tol     !! desired length of the interval of uncertainty
    !! of the final result (>=0)
    real(dl), intent(out)             :: xzero   !! abscissa approximating a zero of f in the interval ax, bx
    real(dl), intent(out)             :: fzero   !! value of f at the root (f(xzero))
    integer, intent(out)              :: iflag   !! status flag (-1 = error, 0 = root found)
    real(dl), intent(in), optional     :: fax     !! if f(ax) is already known, it can be input here
    real(dl), intent(in), optional     :: fbx     !! if f(bx) is already known, it can be input here
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
    real(dl), intent(in)           :: ax      !! left endpoint of initial interval
    real(dl), intent(in)           :: bx      !! right endpoint of initial interval
    real(dl), intent(in)           :: tol     !! desired length of the interval of uncertainty
    !! of the final result (>=0)
    real(dl), intent(out)          :: xzero   !! abscissa approximating a zero of f in the interval ax, bx
    real(dl), intent(out)          :: fzero   !! value of f at the root (f(xzero))
    integer, intent(out)           :: iflag   !! status flag (-1 = error, 0 = root found)
    real(dl), intent(in), optional :: fax     !! if f(ax) is already known, it can be input here
    real(dl), intent(in), optional :: fbx     !! if f(bx) is already known, it can be input here
""",
        )

    def test_alignment_keeps_a_compressible_subblock_before_unaligned_declarations(self) -> None:
        source = """\
    real(dl), intent(in)              :: ax
    real(dl), intent(in)              :: bx
    real(dl), intent(in), optional     :: fax
    real(dl), parameter :: one = 1._dl
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\
    real(dl), intent(in)           :: ax
    real(dl), intent(in)           :: bx
    real(dl), intent(in), optional :: fax
    real(dl), parameter :: one = 1._dl
""",
        )

    def test_declaration_alignment_never_adds_padding(self) -> None:
        source = """\\
    type(c_ptr) :: cptr
    type(CAMBParams), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    class(CAMBParams), target :: this
    type(CAMBParams), pointer :: p
    class(TPythonInterfacedClass) :: replace_with
"""
        self.assertEqual(
            format_text(source, wrap=False),
            """\\
    type(c_ptr) :: cptr
    type(CAMBParams), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    class(CAMBParams), target :: this
    type(CAMBParams), pointer :: p
    class(TPythonInterfacedClass) :: replace_with
""",
        )


class RegressionFixTests(unittest.TestCase):
    def test_leading_continuation_markers_align_text_at_marker(self) -> None:
        source = """\
    if (.not. F%ReadItem(count) .or. count/=this%Count) &
        & call this%Error('TObjectList_LoadState count mismatch (objects must exist before load)')

    generic :: Add => AddString, TNameValueList_AddDouble, &
        &                 TNameValueList_AddReal, TNameValueList_AddInt,&
        &                 TNameValueList_AddLogical
"""
        expected = """\
    if (.not. F%ReadItem(count) .or. count /= this%Count) &
        call this%Error('TObjectList_LoadState count mismatch (objects must exist before load)')

    generic :: Add => AddString, TNameValueList_AddDouble, &
        TNameValueList_AddReal, TNameValueList_AddInt, &
        TNameValueList_AddLogical
"""
        formatted = format_text(source, wrap=False)
        self.assertEqual(formatted, expected)
        self.assertEqual(format_text(formatted, wrap=False), formatted)

    def test_procedure_modifiers_preserve_local_scope_case_handling(self) -> None:
        cases = (
            ("pure function", "function", "result(Y)"),
            ("elemental function", "function", "result(Y)"),
            ("recursive function", "function", "result(Y)"),
            ("pure subroutine", "subroutine", ""),
            ("elemental subroutine", "subroutine", ""),
            ("recursive subroutine", "subroutine", ""),
        )
        for modifier, kind, result_clause in cases:
            with self.subTest(modifier=modifier):
                header = f"{modifier} F(X) {result_clause}".rstrip()
                if kind == "function":
                    source = f"{header}\n    real :: x, y\n    Y = X\nend function F\n"
                    expected_body = "    real :: x, y\n    y = x\n"
                else:
                    source = f"{header}\n    real :: x\n    X = 1\nend subroutine F\n"
                    expected_body = "    real :: x\n    x = 1\n"

                procedures = extract_procedure_cases(source)
                self.assertEqual(len(procedures), 1)
                self.assertEqual(procedures[0].local_names, frozenset({"x", "y"} if kind == "function" else {"x"}))
                formatted = format_text(source, wrap=False)
                self.assertIn(f"{modifier} F(x)" + (" result(y)" if kind == "function" else ""), formatted)
                self.assertIn(expected_body, formatted)
                self.assertIn(f"end {kind} F\n", formatted)

    def test_format_slash_edit_descriptors_are_not_array_constructors(self) -> None:
        source = "10 FORMAT(/)\n20 FORMAT(2X/)\nx=(/1,2/)\n"
        self.assertEqual(
            format_text(source, wrap=False),
            "10 format(/)\n20 format(2X/)\nx = [1, 2]\n",
        )

    def test_cpp_continuation_body_is_never_treated_as_fortran(self) -> None:
        source = "#define SUB subroutine fake(); \\\nend subroutine fake\nprogram p\nend program p\n"
        self.assertEqual(
            [(item.kind, item.name) for item in extract_declared_names(source)],
            [],
        )
        formatted = format_text(source, wrap=False)
        self.assertIn("#define SUB subroutine fake(); \\\nend subroutine fake\n", formatted)
        self.assertNotIn("#define SUB subroutine fake(); \\\n\nend subroutine fake", formatted)

    def test_comment_between_continued_lines_stays_in_one_statement(self) -> None:
        source = """\
subroutine foo(a, &
! explanation
& b)
integer :: a, b
end subroutine foo
"""
        statements = list(_iter_code_statements_with_lines(source))
        self.assertEqual(statements[0].text, "subroutine foo(a, b)")
        self.assertEqual((statements[0].start_line, statements[0].end_line), (0, 2))
        self.assertEqual(extract_procedure_cases(source)[0].local_names, frozenset({"a", "b"}))

    def test_external_macro_argument_sets_exact_case(self) -> None:
        source = "program p\nimplicit none\ninteger :: x\nx=size\nprint *, size\nend program p\n"
        formatted = format_text(source, wrap=False, macro_cases=["SIZE"])
        self.assertIn("x = SIZE", formatted)
        self.assertIn("print *, SIZE", formatted)
        with patch("sys.argv", ["standardize_fortran.py", "-DSIZE=4", "input.f90"]):
            self.assertEqual(parse_args().macro_names, ["SIZE"])

    def test_source_defined_macros_do_not_force_unmatched_case(self) -> None:
        source = "#define SIZE 4\nprint *, size\n"
        self.assertEqual(format_text(source, wrap=False), "#define SIZE 4\nprint *, size\n")

    def test_scope_ranges_refresh_after_terminal_return_removal(self) -> None:
        arguments = ", ".join(["SIZE", *[f"arg{i:02d}" for i in range(30)]])
        source = (
            "subroutine first\n"
            "    return\n"
            "end subroutine first\n"
            f"subroutine second({arguments})\n"
            "    integer :: SIZE\n"
            "    SIZE = 1\n"
            "end subroutine second\n"
        )
        formatted = format_text(source)
        self.assertIn("subroutine second(SIZE,", formatted)
        self.assertNotIn("subroutine second(size,", formatted)

    def test_named_end_reduces_program_unit_depth(self) -> None:
        lines = ["subroutine a\n", "end subroutine a\n", "\n", "\n", "\n", "x=1\n"]
        normalized = formatter.normalize_program_unit_spacing(lines)
        self.assertEqual(normalized[-4:], ["\n", "\n", "\n", "x=1\n"])

    def test_local_type_components_after_module_contains_do_not_leak(self) -> None:
        source = """\
module m
contains
subroutine s
type :: Local
integer :: WeirdCase
end type Local
end subroutine s
end module m
"""
        self.assertEqual(extract_module_variable_names(source), [])
        self.assertEqual(extract_module_variable_types(source), {})

    def test_mapping_defaults_are_real_mappings(self) -> None:
        cases = formatter.FileDeclarationCases({}, {})
        procedure = formatter.ProcedureDeclarationCases(0, 0, {})
        self.assertIsNone(cases.type_procedure_cases.get(("missing", "member")))
        self.assertIsNone(cases.type_component_cases.get(("t", "x")))
        self.assertIsNone(procedure.local_types.get("missing"))

    def test_atomic_write_preserves_symlink(self) -> None:
        with TemporaryDirectory() as directory:
            directory_path = Path(directory)
            target = directory_path / "real.f90"
            link = directory_path / "link.f90"
            target.write_text("old\n")
            link.symlink_to(target.name)
            formatter.write_text_atomic(link, "new\n")
            self.assertTrue(link.is_symlink())
            self.assertEqual(target.read_text(), "new\n")

    def test_external_diff_display_path_does_not_raise(self) -> None:
        path = Path("/tmp/external_source.f90")
        self.assertEqual(formatter._display_path(path), path)

    def test_invalid_extension_is_rejected_before_reading(self) -> None:
        with self.assertRaisesRegex(ValueError, "Expected a free-form Fortran source"):
            formatter._validated_fortran_path(Path("does-not-exist.txt"))

    def test_standard_free_form_extensions_are_accepted(self) -> None:
        for extension in (".f90", ".f95", ".f03", ".f08", ".f18", ".f23", ".F90", ".F23"):
            with self.subTest(extension=extension):
                self.assertEqual(formatter._validated_fortran_path(Path(f"source{extension}")).suffix, extension)

    def test_wrapped_statement_keeps_not_against_its_bracket(self) -> None:
        source = (
            "    subroutine s\n"
            "                if (.not.(State%closed .and. "
            "nint(CTrans%q%points(q_ix)*State%curvature_radius) <= CTrans%ls%l(j))) then\n"
            "                    dlnk = 1\n"
            "                end if\n"
            "    end subroutine s\n"
        )
        formatted = format_text(source)
        self.assertIn("if (.not. (State%closed .and. &\n", formatted)
        self.assertEqual(format_text(formatted), formatted)


if __name__ == "__main__":
    unittest.main()
