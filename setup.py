from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("fastqio.fastq_cython", ["fastqio/fastq_cython.pyx"],
              include_dirs=[np.get_include()])
]

setup(
    name="fastqio",
    version="0.1.0",
    description="FASTQ parser with multi-threaded I/O and Cython-accelerated string handling",
    author="Joseph Chao-Chung Kuo",
    author_email="chaochung.kuo@rwth-aachen.de",
    url="https://github.com/chaochungkuo/fastqio",
    packages=find_packages(),
    ext_modules=cythonize(extensions, language_level=3),
    install_requires=[
        "Cython>=0.29",
        "numpy",
        "pyarrow",
    ],
    python_requires=">=3.7",
)