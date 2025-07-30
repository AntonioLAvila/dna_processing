from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "dna_processing",
        ["dna_processing/dna_processing.cpp"],
        extra_compile_args=["-std=c++20"],
    )
]

setup(
    name="dna_processing",
    version="0.1",
    author="Antonio Avila",
    packages=["dna_processing"],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    include_package_data=True,
    package_data={"dna_processing": ["*.pyi", "py.typed"]},
    zip_safe=False,
)
