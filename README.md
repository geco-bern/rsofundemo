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

To install the development releases of the package run the following commands:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("geco-bern/rsofun")
library("rsofun")
```

### Install manually

If none of the previously mention installation works, download the tar.gz folder from [here](https://github.com/geco-bern/rsofun/releases/tag/v4.4) and use the command

``` r
install.packages(path_to_file, repos = NULL, type="source")
library("rsofun")
```

**NOTE:** Installing from GitHub requires compilation of Fortran and C source code contained in {rsofun}. To enable compiling source code, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or [Xcode](https://developer.apple.com/xcode/) and the [GNU Fortran compiler on Mac](https://github.com/fxcoudert/gfortran-for-macOS) (see also 'Mandatory tools' [here](https://mac.r-project.org/tools/)). On Linux, the gfortran compiler is usually installed already.

Vignettes are not rendered by default, if you want to include additional documentation please use:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("geco-bern/rsofun", build_vignettes = TRUE)
library("rsofun")
```

## Analysis usage

The [00_validation_files_creation.R](https://github.com/stineb/rsofundemo/tree/main/analysis) file creates the validation data necessary to run the P-model. The driver data can be obatined from [FluxDatakit](https://github.com/geco-bern/FluxDataKit/tree/main/data). the files will be stoted in [data](https://github.com/stineb/rsofundemo/tree/main/data).

The file will use two custom function that can be found in [R](https://github.com/stineb/rsofundemo/tree/main/R). To work it requires the path for each file needed and an output path to store the results.

to run the script, you need to use the repository [FLuxDataKit](https://github.com/geco-bern/FluxDataKit/tree/main) to create the files that are needed. The metadata file is already present in [data-raw](https://github.com/stineb/rsofundemo/tree/main/data-raw).

## P-model Use

Below sections show the ease of use of the package in terms of model parameter specification and running both a single run or optimizing the parameters for a given site. For an in depth discussion we refer to the [vignettes](https://geco-bern.github.io/rsofun/articles/).

### Running model

With all data prepared we can run the P-model using `runread_pmodel_f()`. This function takes the nested data structure and runs the model for each site separately, returning nested model output results matching the input drivers. You can use either the `p_model_drivers` and `p_model_validation` which are present in the library or use the data obtained with [00_validation_files_creation.R](https://github.com/stineb/rsofundemo/tree/main/analysis).

the data can be loaded using readRDS, alternatively you can use `p_model_drivers` and `p_model_validation` from `rsofun` package. The rsofun_validation_data.rds file does not contain all the sites (up to now contains only these two)

```{r}
sites <- c("ES-Amo","FR-Pue")
driver <- readRDS(paste0(here("data//"),"rsofun_driver_data.rds"))
validation <- readRDS(paste0(here("data//"),"rsofun_validation_data.rds"))

# filtering (not needed if you run all the sites)
driver <- driver[driver$sitename %in% sites,]
validation <- validation[validation$sitename %in% sites,]
```

## Running model

With all data prepared we can run the P-model using `runread_pmodel_f()`. This function takes the nested data structure and runs the model site by site, returning nested model output results matching the input drivers.

```{r}

# define model parameter values from previous
# work
params_modl <- list(
  kphio              = 0.04998,    # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a        = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b        = 1.0,
  soilm_thetastar    = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao        = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax        = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim         = 30.0,
  kc_jmax            = 0.41
)

# run the model for these parameters
output <- rsofun::runread_pmodel_f(
  driver,
  par = params_modl
)
```

## Plotting creation

TO show the result, a data frame containing all the sites is created.

```{r}
# Create data.frame for plotting
df_gpp_plot <- rbind(
  output |>
    unnest(data) |>
    select(date, gpp, sitename) |>
    mutate(type = "P-model output"),
  driver |>
    unnest(forcing) |>
    select(date, gpp, sitename) |>
    mutate(type = "Observed")
)
df_gpp_plot$type <- factor(df_gpp_plot$type,
                           levels = c('P-model output',
                                      'Observed'))

```

## Output creation

The output is not directly visualized here but is located in [fig](https://github.com/stineb/rsofundemo/tree/main/fig). The model result is shown separately for each site.

```{r}
# Plot GPP

for(i in sites){
  
  p = ggplot(data = df_gpp_plot |> filter(sitename==i)) +
  geom_line(
    aes(x = date,
        y = gpp,
        color = type),
    alpha = 0.7
  ) +
    ggtitle(i)+
  scale_color_manual(values = c(
    'P-model output'='grey70',
    'Observed'='black')) +
  theme_classic() +
  theme(panel.grid.major.y = element_line()) +
  labs(
    x = 'Date',
    y = expression(paste("GPP (g C m"^-2, "s"^-1, ")")),
    colour = ""
  )
  ggsave(paste0(here("fig//"),i, "_p_model_not_calibrated.png"))
}
```

## Calibrating model parameters

To optimize new parameters based upon driver data and a validation dataset we must first specify an optimization strategy and settings, as well as a cost function and parameter ranges.

```{r}
# calibrating
settings <- list(
  method              = "GenSA",
  metric              = cost_rmse_pmodel,
  control = list(
    maxit = 100),
  par = list(
    kphio = list(lower=0.02, upper=0.2, init = 0.05)
  )
)
```

`rsofun` supports both optimization using the `GenSA` and `BayesianTools` packages. The above statement provides settings for a `GenSA` optimization approach. For this example the maximum number of iterations is kept artificially low. In a real scenario you will have to increase this value orders of magnitude. Keep in mind that optimization routines rely on a cost function, which, depending on its structure influences parameter selection. A limited set of cost functions is provided but the model structure is transparent and custom cost functions can be easily written. More details can be found in the "Parameter calibration and cost functions" vignette.

In addition starting values and ranges are provided for the free parameters in the model. Free parameters include: parameters for the quantum yield efficiency `kphio`, `kphio_par_a` and `kphio_par_b`, soil moisture stress parameters `soilm_thetastar` and `soilm_betao`, and also `beta_unitcostratio`, `rd_to_vcmax`, `tau_acclim` and `kc_jmax` (see `?runread_pmodel_f`). Be mindful that with newer versions of `rsofun` additional parameters might be introduced, so re-check vignettes and function documentation when updating existing code.

With all settings defined the optimization function `calib_sofun()` can be called with driver data and observations specified. Extra arguments for the cost function (like what variable should be used as target to compute the root mean squared error (RMSE) and previous values for the parameters that aren't calibrated, which are needed to run the P-model).

```{r}
# calibrate the model and optimize free parameters
pars <- calib_sofun(
  drivers = driver,  
  obs = validation,
  settings = settings,
  # extra arguments passed to the cost function:
  targets = "gpp",             # define target variable GPP
  par_fixed = params_modl[-1]  # fix non-calibrated parameters to previous 
  # values, removing kphio
)
```

The updated model will yield better result.

```{r}
# Update the parameter list with calibrated value
params_modl$kphio <- pars$par["kphio"]

# Run the model for these parameters
output_new <- rsofun::runread_pmodel_f(
  driver,
  par = params_modl
)
```

## Updating the plot dataframe

```{r}
# Update data.frame for plotting
df_gpp_plot <- rbind(
  df_gpp_plot,
  output_new |>
    unnest(data) |>
    select(date, gpp, sitename) |>
    mutate(type = "P-model output (calibrated)")
)
df_gpp_plot$type <- factor(df_gpp_plot$type,
                           levels = c('P-model output',
                                      'P-model output (calibrated)',
                                      'Observed'))
```

## PLotting results

```{r}
for(i in sites){
  
  p = ggplot(data = df_gpp_plot |> filter(sitename==i)) +
  geom_line(
    aes(x = date,
        y = gpp,
        color = type),
    alpha = 0.7
  ) +
  scale_color_manual(values = c(
    'P-model output'='grey70',
    'P-model output (calibrated)'='grey40',
    'Observed'='black')) +
  theme_classic() +
  theme(panel.grid.major.y = element_line()) +
  labs(
    x = 'Date',
    y = expression(paste("GPP (g C m"^-2, "s"^-1, ")")),
    colour = ""
  )
  ggsave(paste0(here("fig//"),i, "_p_model_calibrated.png"))
}
```

## References

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T., and Prentice, I. C.: P-model v1.0: an optimality-based light use efficiency model for simulating ecosystem gross primary production, Geosci. Model Dev., 13, 1545–1581, <https://doi.org/10.5194/gmd-13-1545-2020>, 2020.

For details on the optimization settings we refer to the manuals of [GenSA](https://cran.r-project.org/package=GenSA) and [BayesianTools](https://github.com/florianhartig/BayesianTools).
