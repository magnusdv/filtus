#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'filtus',
    'version': '1.0.1',
    'description': 'Analysis of exome variant files',
    'long_description': open('README.md').read(),
    'author': 'Magnus Dehli Vigeland',
    'author_email': 'magnusdv@medisin.uio.no',
    'license': 'GPL-2',
    'url': 'https://github.com/magnusdv/filtus',
    #'download_url': 'https://github.com/magnusdv/filtus/tarball/v1.0.0',
    'install_requires': [
        'matplotlib',
        'pmw'
    ],
    'packages': [
        'filtus'
    ],
    'classifiers': [
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)'
    ],
    'keywords': [
        'exome',
        'variant filtering',
        'vcf',
        'de novo',
        'autozygosity mapping'
    ],
    'entry_points': {
        'console_scripts': [
            'filtus = filtus.Filtus:main'
        ]
    },
    'package_data': {
        'filtus': [
            'man/pictures/*.png',
            'man/*.html',
            'data/*',
            'testfiles/*'
        ]
    }
}

setup(**config)
