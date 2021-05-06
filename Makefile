.PHONY : install test check build docs clean push_release

test:
	rm -f .coverage
	nosetests --verbose --with-coverage --cover-package ngs_utils tests/*

check:
	flake8 ngs_utils && echo OK
	yapf -r --diff ngs_utils && echo OK

build:
	python setup.py sdist bdist_wheel

docs:
	sphinx-build -a docs docs/_build

clean:
	rm -rf build
	rm -rf dist
	rm -rf ngs_utils.egg-info
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
