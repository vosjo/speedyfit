import setuptools

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
]

setuptools.setup(
    name="speedyfit", # Replace with your own username
    version="0.0.1",
    author="Joris Vos",
    author_email="joris.vos@uv.cl",
    description="MC approach to fit photometric SEDs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vosjo/speedyfit",
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    include_package_data=True,
    test_suite='pytest.collector',
    tests_require=['pytest'],
    entry_points={
        'console_scripts': ['speedyfit=speedyfit.speedyfit:main'],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires='>=3.7',
)