# CSGLIB

This package includes standalone python tools/pipelines for various analyses as well as standalone python modules that may be re-used in other projects. Users are also very welcomed to copy and re-use (cherry-pick) specific code snippets from this package in their personal projects.

For detailed description of available tools and API visit <https://statgen.github.io/csglib/>.

## Installing

### pipenv

pipenv would be the recommended way to install this package if you are developing an application.

You can install directly from github:

```
pipenv install -e git+https://github.com/statgen/csglib.git#egg=csg
```

pipenv specifically recommends installing version-controlled dependencies in editable mode (`-e`).

You can also install a particular commit directly:

```
pipenv install -e git+https://github.com/statgen/csglib.git@08c92d6d846fe364e222d1ea2466f6faac866786#egg=csg
```

### pip

Works basically the same as pipenv above. Note you will need a recent version of pip for this to work (>10). You can upgrade your version of pip by doing `pip install --upgrade pip`.

To install from github:

```
pip install -e git+https://github.com/statgen/csglib.git#egg=csg
```

## Development

Clone the repository:

```
# Clone via ssh
git clone git@github.com:statgen/csglib.git

# Alternatively clone via HTTPS
git clone https://github.com/statgen/csglib.git
```

Create a virtualenv first:

```
virtualenv venv
source venv/bin/activate
```

Then install the package in editable mode with either pip or pipenv:

```
# with pip
pip install -e csglib

# with pipenv
pipenv install -e csglib
```

Pipenv will pick up the surrounding virtualenv if you activate it, rather than creating one in your home directory.

## Building documentation

Updating documentation at gh-pages:

1. On master branch

```
cd docs
sphinx-apidoc -f -o source/ ../csg/
sphinx-build -b html source/ build/
cd ..
git add docs/
git commit -m "documentation updated"
git push
```

2. Modify gh-pages branch

```
git checkout gh-pages
git merge --no-ff -s recursive -X subtree=docs/build/ master
git push
git checkout master
```
