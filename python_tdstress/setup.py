"""
Setup script for python_tdstress package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="python_tdstress",
    version="1.0.0",
    author="Mehdi Nikkhoo (original MATLAB), Python translation team",
    author_email="mehdi.nikkhoo@gmail.com",
    description="Python implementation of triangular dislocation stress calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-repo/TriBIE",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
    ],
)
