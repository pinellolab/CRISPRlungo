#!/usr/bin/env python
from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
        'pandas',
        'matplotlib',
        'seaborn',
        'jinja2',
        'scipy',
        'numpy',
        ]

test_requirements = [ 
        'pandas',
        'matplotlib',
        'seaborn',
        'jinja2',
        'scipy',
        'numpy',
        'unittest',
        ]

__version__ = None
with open('src/CRISPRlungo/CRISPRlungoCore.py') as fin:
    for line in fin:
        if line.startswith('__version__'):
            version_line = line.replace("  "," ")
            version_els = version_line.split(" ")
            __version__ = version_els[2].replace('"','').replace("'","")


setup(
    name='CRISPRlungo',
    author="Kendell Clement",
    author_email='k.clement.dev@gmail.com',
    description="CRISPRlungo enables analysis of single-anchor amplicon sequencing to quantify complex genome editing outcomes.",
    version=__version__,
    entry_points={
        'console_scripts': [
            'CRISPRlungo=CRISPRlungo.cli:crisprlungo',
            'CRISPRlungoBatch=CRISPRlungo.cli:batch',
            'CRISPRlungoCompare=CRISPRlungo.cli:compare',
        ],
    },
    install_requires=requirements,
    long_description=readme,
    include_package_data=True,
    keywords='CRISPRlungo',
    packages=find_packages(where='src',include=['CRISPRlungo', 'CRISPRlungo.*']),
    package_dir={'':'src'},
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/kclem/crisprlungo',
    zip_safe=False,
)
