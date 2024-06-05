from setuptools import setup, find_packages

setup(
    name='QED',
    version='0.1',
    author='JS Seo',
    description='Querying several gene sets to EnrichR at once',
    package_data={
        'qed.data': ['DB/*.csv'],
    },
    install_requires=[
        'tqdm',
        'pandas',
        'requests',
        'scipy',
        'dataclasses;python_version<"3.7"'],
    packages=find_packages(),
)
