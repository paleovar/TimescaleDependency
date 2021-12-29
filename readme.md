# Readme

This repository provides code to reproduce the figures of **"Probing the timescale dependency of local and global variations in surface air temperature from climate simulations and reconstructions of the last millennia"** (Ellerhoff and Rehfeld, 2021) published in *Physical Review E* (2021).

**URL:** https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.064136

**DOI:** 10.1103/PhysRevE.104.064136

**Authors:** Beatrice Ellerhoff and Kira Rehfeld 

**Contact:** beatrice.ellerhoff(at)iup.uni-heidelberg.de

Please see `./license.md` for terms of use. This repository contains the **maintained code and pre-processed data** to create the figures of *Ellerhoff and Rehfeld (2021)*. The raw data will be soon (after the holiday season) archived on [Zenodo](https://www.zenodo.org/), together with the code to generate the pre-processed data. 

# Organisation of this repository

The repository contains scripts (e.g., `F2.R`), input data (`data`), auxillarly files (`helpers`), and information (e.g., `license.md`). To reproduce a figure from *Ellerhoff and Rehfeld (2021)* run the script that is named as the figure (containing its number). Figure `X` of the main manuscript is denoted by `FX.R` and Figure `X` of the supplementary materials is denoted by `FSX.R`. In some cases (e.g. `F3_and_FS11.R`) it was convenient to create several figures together as they build on similar input data. However, sections in the code are marked with *#FX*, so that the individual figures can be also be created from these scripts.  

# Summary

scripts | description
---- | ----------
`F2.R`- `F6.R` | Scripts to reproduce the figures of the **main manuscript**.
`FS1.R`- `FS16.R`| Scripts to reproduce the figures of the **supplementary material**.

directories | description
---- | ----------
`data` | Contains the pre-processed data that serves as input for all figures. The sub-directory `./data/supp` contains data used for supplementary figures only. `./data/shapes/50m_physical` provides mapping information such as coastlines etc. from [naturalearthdata](https://www.naturalearthdata.com/downloads/110m-physical-vectors/), see `ne_50m_land.README.html`. 
`helpers`| Contains scripts (`.R`-files) that define useful functions, initial parameters and load required packages, and meta data (`.Rds` files).

additional files | description
---- | ----------
`.gitignore` | Information for GIT version control to not add several file extensions to version control (e.g. `*.png`, `*.pdf`)
`license.md`/ `license.html` | Licensing information
`readme.md` | General README

# Prerequisites

Our code requires the following [R](https://www.r-project.org/) packages (loaded in `./helpers/init.R`):

- `dplyr`
- `tibble`
- `tidyr`
- `zoo`
- `RColorBrewer`
- `ggplot2`
- `latex2exp`
- `tseries`
- `purrr`
- `PaleoSpec` (can be obtained from `https://github.com/EarthSystemDiagnostics/paleospec`, using `devtools::install_github()`)

Creating some of the figures additionally requires (separately loaded in the corresponding `.R` scripts):

- `cowplot`, `forcats`, `ggnewscale`, `sp`, `viridisLite`, `raster`, `ggcorrplot`, `ggridges`, `grid`, `gridExtra`
- `nest` (can be obtained from `https://github.com/krehfeld/nest`, using `devtools::install_github()`)
- `rgdal` (in case you encounter difficulties installing the package, this might help: https://gist.github.com/dncgst/111b74066eaea87c92cdc5211949cd1e)
- `sf` (in case you encounter difficulties installing the package, this might help: https://r-spatial.github.io/sf/)

**We acknowledge the [R Core team](https://www.R-project.org/) and all package developers of packages used in this study. We thank them for their time and dedication to provide R and the packages to the public. Please see `citation()` for details on the R Core Team and `citation("pkgname")` for details on the developers of individual packages.**


----

The Authors, Dezember, 2021