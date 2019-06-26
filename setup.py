# Lib
from setuptools import setup, find_packages

setup(
    name='methpype',
    version='1.0.0',
    description='Python-based Illumina methylation array preprocessing software',
    long_description='Python-based Illumina methylation array preprocessing software',
    url='https://github.com/LifeEGX/methpype',
    license='MIT',
    author='Life Epigenetics',
    author_email='info@lifeegx.com',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'statsmodels',
        'tqdm'
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    entry_points='''
        [console_scripts]
        methpype-cli=methpype.cli:app
    ''',
)
