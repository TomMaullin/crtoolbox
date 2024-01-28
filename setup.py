from setuptools import setup, find_packages

# Read in the requirements.txt file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="crtoolbox",
    version="0.1.7",
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
    install_requires=requirements, 
    python_requires='>=3.6',
)
