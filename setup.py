from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    "numpy",
    "emcee",
    "astropy",
    "astroquery",
    "pyvo",
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
    version="0.2.0",
    author="Joris Vos",
    author_email="joris.vos@uv.cl",
    description="MC approach to fit photometric SEDs",
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
