#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 17:18:26 2024

@author: tashatu
"""

from setuptools import setup, find_packages
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='SIN_construction',
    version='0.1',
    packages=find_packages(),
    install_requires=requirements,
    author='Tasha Tu',
    author_email='tashatu.bt11@nycu.edu.tw',
    description='A python package to construct single-sample network to analyze individual characteristics.',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url='',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    test_suite='tests',
)