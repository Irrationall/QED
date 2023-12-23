from setuptools import setup, find_packages

setup(
    name='QED',
    version='0.1',
    author='JS Seo',
    description='Querying several gene sets to EnrichR at once',
    install_requires=[
        'tqdm',
        'pandas',
        'requests',
        'dataclasses;python_version<"3.7"'],
    packages=find_packages(),
)