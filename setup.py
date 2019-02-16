#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: maybelinot
# @Email: edik.trott@yandex.ru
# @Date:   2015-09-12 14:07:48
# @Last Modified by:   Eduard Trott
# @Last Modified time: 2015-09-12 21:22:40

from __future__ import absolute_import  # , unicode_literals

from setuptools import setup

VERSION_FILE = "findltr/_version.py"
VERSION_EXEC = ''.join(open(VERSION_FILE).readlines())
__version__ = ''
exec(VERSION_EXEC)  # update __version__
if not __version__:
    raise RuntimeError("Unable to find version string in %s." % VERSION_FILE)

# acceptable version schema: major.minor[.patch][-sub[ab]]
__pkg__ = 'findltr'
__pkgdir__ = {'findltr': 'findltr'}
__pkgs__ = ['findltr']
__provides__ = ['findltr']
__desc__ = 'Get membership details from various data sources.'
__scripts__ = ['bin/findltr', 'bin/ltrfamilies']
__irequires__ = [
    # CORE DEPENDENCIES
    # 'functioncache==0.92', ???
    'bcbio-gff',
    'argparse',
    'pyyaml>=4.2b1',
    'biopython==1.64',
]
__xrequires__ = {
    'tests': [
        'pytest==2.7.2',
        'instructions',
        # 'pytest-pep8==1.0.6',  # run with `py.test --pep8 ...`
    ],
    # 'docs': ['sphinx==1.3.1', ],
    # 'github': ['PyGithub==1.25.2', ],
    # 'invoke': ['invoke==0.10.1', ],

}

pip_src = 'https://pypi.python.org/packages/src'
__deplinks__ = []

# README is in the parent directory
readme_pth = 'README.rst'
with open(readme_pth) as _file:
    readme = _file.read()

github = 'https://github.com/maybelinot/findltr'
download_url = '%s/archive/master.zip' % github

default_setup = dict(
    url=github,
    license='GPLv2',
    author='Eduard Trott',
    author_email='etrott@redhat.com',
    download_url=download_url,
    long_description=readme,
    data_files=[],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Office/Business',
        'Topic :: Utilities',
    ],
    keywords=['information'],
    dependency_links=__deplinks__,
    description=__desc__,
    install_requires=__irequires__,
    extras_require=__xrequires__,
    name=__pkg__,
    package_dir=__pkgdir__,
    packages=__pkgs__,
    provides=__provides__,
    scripts=__scripts__,
    version=__version__,
    zip_safe=False,  # we reference __file__; see [1]
)

setup(**default_setup)
