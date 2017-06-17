# CSGLIB

This package includes standalone python tools/pipelines for various analyses as well as standalone python modules that may be re-used in other projects. Users are also very welcomed to copy and re-use (cherry-pick) specific code snippets from this package in their personal projects.

# Documentation

For detailed description of available tools and API visit <https://statgen.github.io/csglib/>.

# Development

An easy way to install this while developing is: 

```
pip install --user -e .
```

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
