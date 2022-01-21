# Lib
from setuptools import setup, find_packages
exec(open('methylprep/version.py').read())

test_requirements = [
    'methylcheck', # 'git+https://github.com/FoxoTech/methylcheck.git@feature/v0.7.7#egg=methylcheck',
    'pytest',
    'pytest_mock',
    'matplotlib',
    'scikit-learn', # openpyxl uses this, and forcing it to install the best version, not sklearn 0.0
    'openpyxl',
    'coverage'
]

setup(
    name='methylprep',
    version=__version__,
    description='Python-based Illumina methylation array preprocessing software',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    project_urls = {
        "Documentation": "https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/",
        "Source": "https://github.com/FOXOBioScience/methylprep/",
        "Funding": "https://FOXOBioScience.com/"
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Financial and Insurance Industry',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
      ],
    keywords='methylation dna data processing epigenetics illumina',
    url='https://github.com/FOXOBioScience/methylprep',
    license='MIT',
    author='Life Epigenetics',
    author_email='info@FOXOBioScience.com',
    packages=find_packages(),
    include_package_data=True,
    package_data={"":["*.txt.gz"]},
    install_requires=[
        'pyparsing > 3.0',
        'numpy',
        'pandas >=1.3.0',
        'scipy',
        'statsmodels',
        'tqdm',
        'bs4',
        'lxml',
        'requests',
    ],
    extras_require={
        'dev': test_requirements
    },
    setup_requires=['pytest-runner'],
    tests_require= test_requirements,
    entry_points='''
        [console_scripts]
        methylprep-cli=methylprep.cli:app
    ''',
)
