#!/usr/bin/env python3

# Load modules
import os
from setuptools import setup

# Load package info from __version__.py
about = dict()
here = os.path.abspath(os.path.dirname(__file__))
pkg_path = os.path.join(here, "oncopipe")
with open(os.path.join(here, "oncopipe", "__version__.py")) as f:
    exec(f.read(), about)

# Setup package based on package info
setup(
    name=about["__title__"],
    version=about["__version__"],
    description=about["__description__"],
    url=about["__url__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    license=about["__license__"],
    packages=["oncopipe"],
    package_dir={"oncopipe": pkg_path},
    install_requires=["pyyaml", "pandas", "snakemake>=5.31", "packaging"],
    zip_safe=False,
    python_requires=">=3.6.0",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="snakemake bioinformatics workflows",
)
