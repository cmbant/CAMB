# Contributing to CAMB

Thank you for your interest in contributing to CAMB! This guide will help you set up your development environment and understand our code standards.

## Development Setup

### 1. Clone and Install

```bash
git clone --recursive https://github.com/cmbant/CAMB.git
cd CAMB
pip install -e .[dev] [--user]
```

### 2. Install Pre-commit Hooks

```bash
pre-commit install
```

This will automatically format your code and check for issues before each commit.

### 3. Code Formatting Standards

CAMB uses [Ruff](https://docs.astral.sh/ruff/) for code formatting and linting:

- **Line length**: 120 characters
- **Quote style**: Double quotes
- **Python version**: 3.10+
- **Import sorting**: Automatic via ruff

### 4. Before Committing

The pre-commit hooks will run automatically, but you can also run them manually:

```bash
# Run all pre-commit hooks
pre-commit run --all-files
```

## Testing

Run the test suite to ensure your changes don't break anything:

```bash
python -m unittest camb.tests.camb_test
```

## Pull Request Guidelines

1. **Install pre-commit hooks** before making changes
2. **Test your changes** locally
3. **Write clear commit messages**
4. **Keep PRs focused** on a single feature or fix
5. **Update documentation** if needed

## Questions?

- Check the [documentation](https://camb.readthedocs.io/)
- Ask on [CosmoCoffee](https://cosmocoffee.info/viewforum.php?f=11)
- Open an issue for bugs or feature requests

## For AI Agents

If you're an AI agent contributing code:

1. **Always run formatting** before committing:

   ```bash
   ruff format camb fortran/tests
   ruff check --fix camb fortran/tests
   ```

2. **Use clear commit messages** indicating AI authorship:

   ```bash
   git commit -m "Fix issue XYZ (AI agent)"
   ```

3. **Test thoroughly** before submitting PRs. You may need to install gfortran first, and possibly recursively clone the
   forutils submodule.
