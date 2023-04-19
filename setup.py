from setuptools import setup, find_packages
import codecs
import os
import re


def read_meta():
    """
    Read the meta information from the nowcast.__init__.py file.  Assume UTF-8 encoding.
    """
    meta_path = os.path.join("speedyfit", "__init__.py")
    with codecs.open(os.path.join(os.path.abspath(os.path.dirname(__file__)), meta_path), "rb", "utf-8") as f:
        return f.read()


def find_meta(meta):
    """
    Extract __*meta*__ from META_FILE.
    """
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta),
        META_FILE, re.M
    )
    if meta_match:
        return meta_match.group(1)
    raise RuntimeError("Unable to find __{meta}__ string.".format(meta=meta))


META_FILE = read_meta()

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    "numpy",
    "emcee",
    "astropy>=5.1",
    "astroquery",
    "pyvo>=1.4",
    "matplotlib",
    "pyaml",
    "corner",
    "tqdm",
    "scipy",
    "h5py",
    "pandas",
    "gaiadr3-zeropoint",
]
setup(
    name="speedyfit",
    version=find_meta("version"),
    author=find_meta("author"),
    author_email=find_meta("email"),
    description=find_meta("description"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vosjo/speedyfit",
    packages=find_packages(),
    include_package_data=True,
    install_requires=install_requires,
    test_suite='pytest.collector',
    tests_require=['pytest'],
    entry_points={
        'console_scripts': ['speedyfit=speedyfit.main:main', 'speedyfit-batch=speedyfit.speedyfit_batch:main'],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires='>=3.7',
)
