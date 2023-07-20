# The Confidence Regions Toolbox

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
 - [File Handling](#file-handling)
 - [Regression](#regression) 
 - [Generating CRs](#generating-crs)
 - [Working with Boundaries](#working-with-boundaries)
 - [Displaying Volumes](#displaying-volumes)
- [Data Simulation](#data-simulation)
 - [2D Data Simulation](#2d-data-simulation)
 - [3D Data Simulation](#3d-data-simulation)
- [License](#license)

## Background

The crtoolbox is a Python-based software package designed to provide Confidence Regions (CRs) for excursion sets derived from neuroimaging data. At present, the methods provided are those proposed in [Sommerfeld, et al. 2017](https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1341838), [Bowring, et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33166643/) and [Maullin-Sapey, et al. 2022](https://arxiv.org/abs/2201.02743). These methods provide CRs for 2D and 3D regression coefficient images, Cohen's d effect size images, and conjunction or disjunction images. For each approach, estimated inner and outer sets are constructed for an excursion set of interest. These estimated sets, together known as Confidence Regions, are designed to bound the true excursion set with a user-specified level of confidence.

Further detail on CRs, alongside an extensive demonstration of the `crtoolbox`, can be found in the [`crtoolbox_demo`](https://github.com/TomMaullin/crtoolbox_demo/) repository, which originally formed part of the `Beyond Blobology` course at the Organization for Human Brain Mapping conference 2023. **Please Note:** The function calls described in this practical have only been tested with the `crtoolbox==0.1.6` version. Later versions may not be compatible.

## Install

To install the `crtoolbox` please use pip:

```
pip install crtoolbox
```

## Usage

At present the `crtoolbox` consists of several standalone functions, which may be used together to generate CRs. The most relevant of these are detailed below. Full docstrings for each function can be found in the files linked.

### File Handling

The following functions can be used for file manipulation and handling:

 - `read_image`: This function may be used to read 2D or 3D images. Supported formats include; RGB images supported by the `Image` package, `numpy` data files and and NIfTI images.
 - `addBlockToNifti`: This function adds a block of voxels to a pre-existing NIFTI or creates a new NIFTI of specified dimensions if the file does not already exist.
 - `remove_files`: This function removes lists of files.

These functions may be found in [`crtoolbox.lib.fileio`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/lib/fileio.py).

### Effect Size images

The following functions can be used for regression and effect size estimation for either 2d or 3d images.

 - `regression`: This function can be used to perform a standard mass-univariate linear regression and obtain beta coefficient estimates ("raw" effect sizes). Please note that for large analyses, other packages may be more suitable ([for instance...](https://github.com/TomMaullin/BLM)).
 - `cohens`: This function can be used to generate cohen's D images.

Please note that the residuals and standard deviation images output by `regression` and `cohens` are not the same; `regression` outputs regression residuals and the MLE for the standard deviation, whilst `cohens` outputs transformed residuals (generated according to method 2 of [Bowring, et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33166643/)) alongside their standard deviation.

These functions may be found in [`crtoolbox.lib.regression`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/lib/regression.py) and [`crtoolbox.lib.cohens`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/lib/cohens.py), respectively.

### Generating CRs

The following function is used to generate confidence regions. It supports both 2d and 3d images.

 - `generate_CRs`: This function uses the outputs of `regression`, or `cohens`, to generate CRs for raw or standardized effect sizes, respectively.

This function may be found in [`crtoolbox.generate`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/generate.py).

### Coverage Assessment

The below function may be used to check whether the CRs correctly bound the underlying true excursion set (note this is only useful for simulation settings). It uses both a binary inclusion check and the interpolation methods outlines in [Bowring, et al. 2019](https://www.sciencedirect.com/science/article/pii/S1053811919307785).

 - `check_violations`: Given a pair of CRs and a "true" excursion set, this function assesses whether the CRs correctly enclose the excursion set.

This function may be found in [`crtoolbox.coverage`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/coverage.py).

### Working with Boundaries

The below function may be useful for developing the tool.

 - `get_bdry_map_combined`: This function produces an image of the boundary of an excursion set. Both 2D and 3D images are supported.

This function may be found in [`crtoolbox.lib.boundary`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/lib/boundary.py).

### Displaying Volumes

The below functions can be used to display various images in an interactive format.

display_volume, display_crs

 - `display_volume`: This function can be used to dispay individual 2D or 3D volumes. Masking is supported.
 - `display_crs`: This function can be used to dispay 2D or 3D CRs. Masking is supported.

These functions may be found in [`crtoolbox.lib.display`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/lib/display.py).

## Data Simulation

The `crtoolbox` currently contains custom code for generating simulated data which may be useful for testing and demonstration. The below sections provide detail on the available 2D and 3D simulated datasets respectively.

### 2D Data Simulation

The below classes can be used to generate 2D signals:

 - `CircleSignal`: This will create a circular signal with user specified radius, center, smoothness and peak magnitude.
 - `SquareSignal`: This will create a square signal with user specified radius, center, smoothness and peak magnitude.
 - `RampSignal`: This will create a linear ramp signal with user specified slope and value range.
 
The above classes will create objects representing each of the signal types. An image of the signal can be viewed using the `plot` method. Noise can also be generated using the `Noise` class.
 
 - `Noise`: This will create noise with user specified magnitude, smoothness and heterogeniety settings.
 
To use the above classes to generate data the `generate_data_2D` function can be used:

 - `generate_data_2D`: This will take in a signal object and a noise object and return the user's specified number of data instances.

The classes for signal generation and noise generation, as well as the `generate_data_2D` function, may be found in [`crtoolbox.tests.generate_2d_data`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/tests/generate_2d_data.py).

### 3D Data Simulation

The below function can be used to generate 3D simulated data resembling real-world fMRI images.

 - `generate_data`: This will generate simualted 3D fMRI images according to a mass-univariate linear regression model, with fixed masking.

For further detail on the pipeline used to generate these images, see [Maullin-Sapey & Nichols. 2022](https://www.biorxiv.org/content/10.1101/2022.03.09.483645v1.full). The `generate_data` function may be found in [`crtoolbox.tests.generate_ni_data`](https://github.com/TomMaullin/crtoolbox/blob/master/crtoolbox/tests/generate_ni_data.py).

## License

[MIT](LICENSE) 
