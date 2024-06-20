"""
scTail: single-cell alternative PAS detection and analysis library
See: https://github.com/StatBiomed/scTail
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path






# Set __version__ for the project.
exec(open("./scTail/version.py").read())

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['torch>=2.1.0','pyranges>=0.0.129','pandas>=2.1.1','anndata>=0.10.3',
'pysam>=0.22.0','scipy>=1.12.0','scikit-learn>=1.3.2','matplotlib>=3.8.1','pyfaidx>=0.7.2.2','tqdm>=4.66.1',
'kipoiseq>=0.7.1','numpy>=1.22.4','attrs<=21.4.0']


setup(
    name='scTail',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='scTail: Detection alternative PAS in single cells',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/StatBiomed/scTail',

    # Author details
    author=['Ruiyan Hou'],
    author_email='ruiyan@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['Polyadenylation sites', 'single-cell RNA-seq'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().




    packages=find_packages(),

    #package_data={'scTail': ['model/*.sav']},

    entry_points={
          'console_scripts': [
            'scTail = scTail.bin.scTail_main:main',
            'scTail-callPeak = scTail.bin.callPeak:main',
            'scTail-peakMerge = scTail.bin.peakMerge:main',
            'scTail-count = scTail.bin.count:main',
            ],
          }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html

    python_requires='>=3.9',

    install_requires=reqs,

    extras_require={
        'docs': [
            #'sphinx == 1.8.3',
            'sphinx_bootstrap_theme']},

    #py_modules = ['scTail']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)
