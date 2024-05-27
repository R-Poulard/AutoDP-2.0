from setuptools import setup

setup(
    name="autoDP",
    packages = ['autodp'],
    package_dir = {'autodp': 'python'},
    py_modules = ['minimal_expansion','tree_decomposition','tree_of_bags']
)
