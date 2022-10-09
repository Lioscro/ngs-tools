from setuptools import find_packages, setup


def read(path):
    with open(path, 'r') as f:
        return f.read()


long_description = read('README.md')

setup(
    name='ngs-tools',
    version='1.8.1',
    url='https://github.com/Lioscro/ngs-tools',
    author='Kyung Hoi (Joseph) Min',
    author_email='phoenixter96@gmail.com',
    description='Reusable tools for working with next-generation sequencing (NGS) data',  # noqa
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='',
    python_requires='>=3.6',
    license='MIT',
    packages=find_packages(exclude=('tests', 'docs', 'tests.*')),
    zip_safe=False,
    include_package_data=True,
    install_requires=read('requirements.txt').strip().split('\n'),
    extras_require={
        'levenshtein': read('requirements-levenshtein.txt').strip().split('\n'),
    },
    entry_points={},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities',
    ],
)
