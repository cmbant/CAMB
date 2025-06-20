# Contributing to CAMB

Thank you for your interest in contributing to CAMB! This guide will help you set up your development environment and understand our code standards.

## Development Setup

### 1. Clone and Install

```bash
git clone --recursive https://github.com/cmbant/CAMB.git
cd CAMB
pip install -e .[dev]
```

### 2. Install Pre-commit Hooks

```bash
pre-commit install
```

This will automatically format your code and check for issues before each commit.

### 3. Code Formatting Standards

CAMB uses [Ruff](https://docs.astral.sh/ruff/) for Python formatting and linting:

- **Line length**: 120 characters
- **Quote style**: Double quotes
- **Python version**: 3.10+
- **Import sorting**: Automatic via ruff
- **Target version**: py310

The pre-commit hooks include:

- **Ruff formatting and linting** (Python files)
- **Trailing whitespace removal** (Python, Fortran, Jupyter notebooks)
- **End-of-file fixing** (Python, Fortran, Jupyter notebooks)
- **PyUpgrade** for Python 3.10+ syntax

### 4. Before Committing

The pre-commit hooks will run automatically, but you can also run them manually:

```bash
# Run all pre-commit hooks
pre-commit run --all-files

# Run only ruff formatting
ruff format camb fortran/tests

# Run only ruff linting with fixes
ruff check --fix camb fortran/tests
```

## Testing

Run the test suite to ensure your changes don't break anything:

```bash
python -m unittest camb.tests.camb_test
```

For HMcode tests (Linux only):

```bash
git clone https://github.com/alexander-mead/HMcode_test_outputs.git
python -m unittest camb.tests.hmcode_test
```

## VS Code Setup

The repository includes VS Code configuration files with:

- **Recommended extensions**: Python, Ruff, Fortran linter
- **Format on save**: Enabled with Ruff as the default Python formatter
- **Rulers at 120 characters**
- **Pylance settings**: Configured to silence most NumPy-related type errors
- **Fortran formatting**: Configured with findent

### Recommended Extensions

The following extensions will be suggested when you open the project:

- `ms-python.python` - Python support
- `ms-python.debugpy` - Python debugging
- `charliermarsh.ruff` - Ruff formatter and linter
- `fortran-lang.linter-gfortran` - Fortran support

### Optional Fortran Tools

For enhanced Fortran development, you may want to install:

```bash
pip install findent fortls
```

Note: You may need to configure global VS Code settings for Fortran tool paths.

## Pull Request Guidelines

1. **Install pre-commit hooks** before making changes
2. **Test your changes** locally
3. **Write clear commit messages**
4. **Keep PRs focused** on a single feature or fix
5. **Update documentation** if needed
6. **Ensure CI passes** - GitHub Actions will check formatting and run tests

## Questions?

- Check the [documentation](https://camb.readthedocs.io/)
- Ask on [CosmoCoffee](https://cosmocoffee.info/viewforum.php?f=11)
- Open an issue for bugs or feature requests

## For AI Agents

If you're an AI agent contributing code, follow these guidelines:

### 1. Use Standard Ruff/Pre-commit Settings

### 2. Clear Attribution in PRs

- **Use author name "AI agent"** when creating PRs
- **Include clear commit messages** indicating AI authorship:
  ```bash
  git commit -m "Fix issue XYZ (AI agent)"
  ```
