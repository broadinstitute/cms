try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

EXCLUDE_FROM_PACKAGES = ['old','cms.test*']

config = {
    'description': 'Composite of Multiple Signals: tests for selection in meiotically recombinant populations',
    'author': 'Broad Institute',
    'url': 'https://github.com/broadinstitute/cms',
    'download_url': 'https://github.com/broadinstitute/cms',
    'author_email': '',
    'version': '0.1',
    'install_requires': ['nose','biopython','pysam'],
    'packages': find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    'scripts': ['cms/selection.py','cms/main.py'],
    'name': 'python-cms'
}

setup(**config)
