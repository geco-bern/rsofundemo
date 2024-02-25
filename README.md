# Demonstration workflows with rsofun

Contains demonstration workflows for using {rsofun} for common applications. This is provided here in a separate package and repository to keep code and dependencies in {rsofun} slim.

All demonstration workflows are implemented as vignettes in `./vignettes/` and displayed under 'Articles' on the [website](https://geco-bern.github.io/rsofun/articles/).

## Installation

### Stable release

**WARNING: rsofun is not currently available on CRAN.** We're working on it. Until it's available again, the command below will not work.

``` r
install.packages("rsofun")
library("rsofun")
```

### Development release

To install the development releases of the package run the following
commands:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("geco-bern/rsofun")
library("rsofun")
```

### Install manually

If none of the previuously mention installation works, download the tar.gz folder from [here](https://github.com/geco-bern/rsofun/releases/tag/v4.4) and use the command

``` r
install.packages(path_to_file, repos = NULL, type="source")
library("rsofun")
```

**NOTE:** Installing from GitHub requires compilation of Fortran and C source code contained in {rsofun}. To enable compiling source code, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or [Xcode](https://developer.apple.com/xcode/) and the [GNU Fortran compiler on Mac](https://github.com/fxcoudert/gfortran-for-macOS) (see also 'Mandatory tools' [here](https://mac.r-project.org/tools/)). On Linux, the gfortran compiler is usually installed already.

Vignettes are not rendered by default, if you want to include additional
documentation please use:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("geco-bern/rsofun", build_vignettes = TRUE)
library("rsofun")
```

## Analysis usage

The [00_validation_files_creation.R](https://github.com/stineb/rsofundemo/tree/main/analysis) file creates the validation data necessary to run the P-model. The driver data can be obatined from [Zenodo](https://zenodo.org/records/8403081).
the files will be stored in [data](https://github.com/stineb/rsofundemo/tree/main/data).

The file will use two custom function that can be found in [R](https://github.com/stineb/rsofundemo/tree/main/R). To work it requires the path for each file needed and an output path to store the results. 

to run the script, you need to download the data from [Zenodo](https://zenodo.org/records/8403081). The metadata file is already present in [data-raw](https://github.com/stineb/rsofundemo/tree/main/data-raw).

## P-model Use

Below sections show the ease of use of the package in terms of model parameter specification and running both a single run or optimizing the parameters for a given site. For an in depth discussion we refer to the [vignettes](https://geco-bern.github.io/rsofun/articles/).

### Running model

With all data prepared we can run the P-model using `runread_pmodel_f()`. This function takes the nested data structure and runs the model for each site separately, returning nested model output results matching the input drivers. The data structure needs to have a specific strucutre explained in `p_model_drivers` documentation.
You can use either the `p_model_drivers` and `p_model_validation` or use the driver data downloaded previously. The validation data can be obtainted from  [00_validation_files_creation.R](https://github.com/stineb/rsofundemo/tree/main/analysis). The last two data can be loaded using `readRDS()`.

## Running model

To run the model, is necessary to define the parameters that will be using during the analysis. The optimal parameters are alreday present in the vignette.
With all data prepared we can run the P-model using `runread_pmodel_f()`. This function takes the nested data structure and runs the model site by site, returning nested model output results matching the input drivers.
To have an in-depth description of the model output refers to `run_pmodel_f_bysite ()` documentation.

## Plotting creation

To show the result, a data frame containing all the sites is created. The plot will be stored as png file which can be found in [fig](https://github.com/stineb/rsofundemo/tree/main/fig). 

## Calibrating model parameters

To optimize new parameters based upon driver data and a validation dataset we must first specify an optimization strategy and settings, as well as a cost function and parameter ranges.

`rsofun` supports both optimization using the `GenSA` and `BayesianTools` packages. The above statement provides settings for a `GenSA` optimization approach. For this example the maximum number of iterations is kept artificially low. In a real scenario you will have to increase this value orders of magnitude. Keep in mind that optimization routines rely on a cost function, which, depending on its structure influences parameter selection. A limited set of cost functions is provided but the model structure is transparent and custom cost functions can be easily written. More details can be found in the "Parameter calibration and cost functions" vignette.

In addition starting values and ranges are provided for the free parameters in the model. Free parameters include: parameters for the quantum yield efficiency `kphio`, `kphio_par_a` and `kphio_par_b`, soil moisture stress parameters `soilm_thetastar` and `soilm_betao`, and also `beta_unitcostratio`, `rd_to_vcmax`, `tau_acclim` and `kc_jmax` (see `?runread_pmodel_f`). Be mindful that with newer versions of `rsofun` additional parameters might be introduced, so re-check vignettes and function documentation when updating existing code.

With all settings defined the optimization function `calib_sofun()` can be called with driver data and observations specified. Extra arguments for the cost function (like what variable should be used as target to compute the root mean squared error (RMSE) and previous values for the parameters that aren't calibrated, which are needed to run the P-model).

The updated model will yield better result. The model will run in the same way as the non calibrated p model.

## Updating the plot dataframe

After the calibration is possible to update the plot dataframe by adding the parameterized output and the plot will be stored as png file.

## References

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T., and Prentice, I. C.: P-model v1.0: an optimality-based light use efficiency model for simulating ecosystem gross primary production, Geosci. Model Dev., 13, 1545–1581, https://doi.org/10.5194/gmd-13-1545-2020, 2020.

For details on the optimization settings we refer to the manuals of [GenSA](https://cran.r-project.org/package=GenSA) and [BayesianTools](https://github.com/florianhartig/BayesianTools).
