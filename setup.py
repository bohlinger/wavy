#!/usr/bin/env python

import setuptools

setuptools.setup(
    name='wavy',
    version='0.0',
    description= 'A package for processing wave measurements and wave model output',
    author='Patrik Bohlinger',
    author_email='patrikb@met.no',
    url='https://github.com/bohlinger/wavy',
    packages=setuptools.find_packages(),
    include_package_data=False,
    setup_requires=['setuptools_scm'],
    )
