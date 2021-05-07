.PHONY : install test check build docs clean push_release

test:
	rm -f .coverage
	nosetests --verbose --with-coverage --cover-package ngs_tools tests/*

check:
	flake8 ngs_tools && echo OK
	yapf -r --diff ngs_tools && echo OK

build:
	python setup.py sdist

docs:
	sphinx-build -a docs docs/_build

clean:
	rm -rf build
	rm -rf dist
	rm -rf ngs_tools.egg-info
	rm -rf docs/_build
	rm -rf docs/api
	rm -rf .coverage

bump_patch:
	bumpversion patch

bump_minor:
	bumpversion minor

bump_major:
	bumpversion major

push_release:
	git push && git push --tags
