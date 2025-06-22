from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="crtoolbox",
    version="0.2.0",
    author="Tom Maullin",
    author_email="TomMaullin@gmail.com",
    description="The Confidence Regions Toolbox",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tommaullin/crtoolbox",
    packages=find_packages(),
    package_data={
        'crtoolbox.tests': ['mask.nii'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "nibabel",
        "nilearn",
        "pyyaml",
        "plotly",
        "nbformat>=4.2.0",
    ],
    python_requires='>=3.6',
)
