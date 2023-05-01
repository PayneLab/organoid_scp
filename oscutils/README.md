# oscutils

Utility package for working with tables output by MetaMorpheus and Proteome Discoverer in the `organoid_scp` project.

## Installation and use

To use this package, install it using `pip`. Use the command `pip install [PATH TO THIS DIRECTORY]`. E.g. if this directory is your current working directory, use 

```
pip install .
```

since `.` is a reference to the current working directory. If you're in the root directory of the repository, use

```
pip install oscutils/
```

Why is this directory special, you may ask? Because this is where the [`pyproject.toml`](https://github.com/PayneLab/organoid_scp/blob/main/oscutils/pyproject.toml) file is, which contains all the metadata that `pip` needs to build the package.

`pip` will automatically read the package metadata to handle software dependencies. However, if you'd like to look at them yourself, they're in [`pyproject.toml`](https://github.com/PayneLab/organoid_scp/blob/main/oscutils/pyproject.toml). At the time this README was last updated, we required Python 3.10 or higher.

After installing the package, you can import and use it in any Python script or notebook using `import oscutils`.

## Style guide when adding code to the package

Functions should be annotated with type hints following the format outlined in [PEP 484](https://peps.python.org/pep-0484/). Note that we use `|` as an abbreviation for `typing.Union`, hence this package requires Python 3.10 or higher.

Docstrings should follow the [`numpydoc` standard](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard).
