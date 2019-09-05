# Lib
from setuptools import setup, find_packages

setup(
    name='methylprep',
    version='1.1.0',
    description='Python-based Illumina methylation array preprocessing software',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    project_urls = {
        "Documentation": "https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/",
        "Source": "https://github.com/lifeEGX/methylprep/",
        "Funding": "https://lifeegx.com/"
    },
    classifiers=[
        'Development Status :: 4 - Beta',
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
    url='https://github.com/LifeEGX/methylprep',
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
        'tqdm',
        'bs4',
        'lxml'
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    entry_points='''
        [console_scripts]
        methylprep-cli=methylprep.cli:app
    ''',
)
