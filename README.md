# Demonstration workflows with rsofun

Contains demonstration workflows for using {rsofun} for common applications. This is provided here in a separate package and repository to keep code and dependencies in {rsofun} slim.

All demonstration workflows are implemented as vignettes in `./vignettes/` and displayed under 'Articles' on the website.


## Installation

### stable release

**WARNING: rsofun is not currently available on CRAN.** We're working on it. Until it's available again, the command below will not work.

``` r
install.packages("rsofun")
library("rsofun")
```

### development release

To install the development releases of the package run the following
commands:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("geco-bern/rsofun")
library("rsofun")
```

### install manually

If none of the previuously mention installation works, download the tar.gz folder from here https://github.com/geco-bern/rsofun/releases/tag/v4.4 and use the command

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

## analysis usage

The [00_p_model_files_creation.R](https://github.com/FrancescoGrossi-unimi/rsofundemo/tree/main/analysis) file creates the data necessary to run the P-model. The data will be stored in rda file format. Two pairs of data (drivers and validation) of two different sites can be found in [files](https://github.com/FrancescoGrossi-unimi/rsofundemo/tree/main/vignettes/files).

If you want to run the P-model with different sites, be sure to have the following files:

- The csv and lsm files of the site downloable [fluxnet.org](https://fluxnet.org/data/fluxnet2015-dataset/)
- the nc and metadata files present in the repo [ancillary data]((https://github.com/FrancescoGrossi-unimi/rsofundemo/tree/main/vignettes/ancillary_data))

NC files too large TO SOLVE

## P-model Use

Below sections show the ease of use of the package in terms of model parameter specification and running both a single run or optimizing the parameters for a given site. For an in depth discussion we refer to the [vignettes](https://geco-bern.github.io/rsofun/articles/).

### Running model

With all data prepared we can run the P-model using `runread_pmodel_f()`. This function takes the nested data structure and runs the model site by site, returning nested model output results matching the input drivers. You can use either the `p_model_drivers` and `p_model_validation` which are present in the library or use the data obtained with [00_p_model_files_creation.R](https://github.com/FrancescoGrossi-unimi/rsofundemo/tree/main/analysis).

```r
# load the data
site <- "FR-Pue"
driver <- readRDS(paste0(path_file,site,"_p_model_drivers.rda"))
validation <- readRDS(paste0(path_file,site,"_p_model_validation.rda"))
```

``` r
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
  p_model_drivers, # or use drivers
  par = params_modl
  )
```

### Parameter optimization

To optimize new parameters based upon driver data and a validation dataset we must first specify an optimization strategy and settings, as well as a cost function and parameter ranges.

``` r
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

``` r
# calibrate the model and optimize free parameters
pars <- calib_sofun(
    drivers = p_model_drivers,  # or use drivers
    obs = p_model_validation, # or use validation
    settings = settings,
    # extra arguments passed to the cost function:
    targets = "gpp",             # define target variable GPP
    par_fixed = params_modl[-1]  # fix non-calibrated parameters to previous 
                                 # values, removing kphio
  )
```

## P-model multisite 

the [rsofun_multisite_fdk.Rmd](https://github.com/FrancescoGrossi-unimi/rsofundemo/tree/main/vignettes) follow the same workflow as P-model. After loading the files is possible to run the model as above

``` r
# output path for images
out_path <- paste0(getwd(),"/vigenttes/img/") # path output for the images
path  <- paste0(getwd(),"/vigenttes/files/")     #instert your path of .rda files 
sites <- c("ES-Amo","FR-Pue")
drivers <-NULL
for(i in sites){          
  tmp <-  readRDS(paste0(path,i,"_p_model_drivers.rda"))
  drivers <- rbind(drivers,tmp)
}
validations <-NULL
for(i in sites){
  tmp <-  readRDS(paste0(path,i,"_p_model_validation.rda"))
  validations <- rbind(validations,tmp)
}
``` 

## parameter initialization and running model 

```r
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
```
```r
# run the model for these parameters
output <- NULL
for (i in sites){
  tmp <- rsofun::runread_pmodel_f(
  drivers |> filter(sitename==i),
  par = params_modl)
  output <- rbind(output,tmp)
 } 
```

## Showing results

The data will not be shown in the vignette, instead a png and PDF files will be created for each site in the folder [img](https://github.com/FrancescoGrossi-unimi/rsofundemo/tree/main/vignettes).

### Data frame creation for plotting

```r
# Create data.frame for plotting
df_gpp_plot <- rbind(
  output |>
    unnest(data) |>
    select(date, gpp, sitename) |>
    mutate(type = "P-model output"),
  validations |>
    unnest(data) |>
    select(date, gpp, sitename) |>
    mutate(type = "Observed")
)
df_gpp_plot$type <- factor(df_gpp_plot$type,
                           levels = c('P-model output',
                                      'Observed'))
```

### Plotting

```r
# Plot GPP
for(i in sites){
  cairo_pdf(filename = paste0(out_path,i, "_p_model_not_calibrated.pdf"), height=5, width=10)
  png(filename = paste0(out_path,i, "_p_model_not_calibrated.png"), height=700*3, width=700*3, res = 300)
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
  p |> print()
  dev.off()
}
```
## Paramter optimization

The paramter are optimized for each site separately follwing the same procedure above

```r
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

```r
# calibrate the model and optimize free parameters

pars <- NULL

for(i in sites){
  tmp <- calib_sofun(
  drivers = drivers |> filter(sitename==i),  
  obs = validations |> filter(sitename==i),
  settings = settings,
  # extra arguments passed to the cost function:
  targets = "gpp",             # define target variable GPP
  par_fixed = params_modl[-1]  # fix non-calibrated parameters to previous 
  # values, removing kphio
  )
  names(tmp) <- paste0(names(tmp),"_",i)
  pars <- append(pars,tmp)
}
```

### Updating

```r
# Update the parameter list with calibrated value

output_new <- NULL
for (i in sites){
  
  params_modl$kphio <- pars[names(pars) == paste0("par_",i)][[1]]
  tmp <- rsofun::runread_pmodel_f(
  drivers|> filter(sitename==i),
  par = params_modl)
  output_new <- rbind(output_new,tmp)
 } 
```

```r
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
### Plotting

```r
# Plot GPP

for(i in sites){
  cairo_pdf(filename = paste0(out_path,i, "_p_model_calibrated.pdf"), height=5, width=10)
  png(filename = paste0(out_path,i, "_p_model_calibrated.png"), height=700*3, width=700*3, res = 300)
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
  p |> print()
  dev.off()
}
```

## References

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T., and Prentice, I. C.: P-model v1.0: an optimality-based light use efficiency model for simulating ecosystem gross primary production, Geosci. Model Dev., 13, 1545–1581, https://doi.org/10.5194/gmd-13-1545-2020, 2020.
