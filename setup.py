from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "dna_processing",
        ["dna_processing.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-std=c++20"]
    ),
]

setup(
    name="dna_processing",
    ext_modules=ext_modules,
)
