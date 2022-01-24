# Readme

This repository provides code to reproduce the figures of **"Probing the timescale dependency of local and global variations in surface air temperature from climate simulations and reconstructions of the last millennia"** (Ellerhoff and Rehfeld, 2021) published in *Physical Review E* (2021).

**URL:** https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.064136

**DOI:** 10.1103/PhysRevE.104.064136

**Authors:** Beatrice Ellerhoff and Kira Rehfeld 

**Contact:** beatrice.ellerhoff(at)iup.uni-heidelberg.de

Please see `./license.md` for terms of use. This repository contains the **maintained code and pre-processed data** to create the figures of *Ellerhoff and Rehfeld (2021)*. The raw data will be soon archived on [Zenodo](https://www.zenodo.org/), together with the code to generate the pre-processed data. It is therefore important to note, that the subfolder (`./processing`) is not complete yet. However, all figure can be reproduced using the post-processed data (`./data`).

## Organisation of this repository

The repository contains scripts (e.g., `F2.R`), input data (`data`), auxillarly files (`helpers`), and information (e.g., `license.md`). To reproduce a figure from *Ellerhoff and Rehfeld (2021)* run the script that is named as the figure (containing its number). Figure `X` of the main manuscript is denoted by `FX.R` and Figure `X` of the supplementary materials is denoted by `FSX.R`. In some cases (e.g. `F3_and_FS11.R`) it was convenient to create several figures together as they build on similar input data. However, sections in the code are marked with *#FX*, so that the individual figures can be also be created from these scripts.  

scripts | description
---- | ----------
`F2.R`- `F6.R` | Scripts to reproduce the figures of the **main manuscript**.
`FS1.R`- `FS16.R`| Scripts to reproduce the figures of the **supplementary material**.

directories | description
---- | ----------
`./data/` | Contains the pre-processed data that serves as input for all figures. The sub-directory `./data/supp` contains data used for supplementary figures only. `./data/shapes/50m_physical` provides mapping information such as coastlines etc. from [naturalearthdata](https://www.naturalearthdata.com/downloads/110m-physical-vectors/), see `ne_50m_land.README.html`. 
`./helpers/`| Contains scripts (`.R`-files) that define useful functions, initial parameters and load required packages, and meta data (`.Rds` files).

additional files | description
---- | ----------
`.gitignore` | Information for GIT version control to not add several file extensions to version control (e.g. `*.png`, `*.pdf`)
`license.md`/ `license.html` | Licensing information
`readme.md` | General README

## Prerequisites

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

- `cowplot`, `forcats`, `ggnewscale`, `sp`, `viridisLite`, `raster`, `ggcorrplot`, `ggridges`, `grid`, `gridExtra`, `irr`
- `nest` (can be obtained from `https://github.com/krehfeld/nest`, using `devtools::install_github()`)
- `rgdal` (in case you encounter difficulties installing the package, this might help: https://gist.github.com/dncgst/111b74066eaea87c92cdc5211949cd1e)
- `sf` (in case you encounter difficulties installing the package, this might help: https://r-spatial.github.io/sf/)

Generating some of the pre-processed data in `./processing` (using data from Zenodo) might further require:
- `ncdf4`
- `foreach`
- `poweRlaw`
- `VGAM`
- `R.matlab`
- `doParallel`
- `furrr`

## Data references

This manuscript is based upon data provided by the World Climate Research Programme’s Working Group on Coupled Modelling, which is responsible for CMIP and
PMIP. We give detailed information on all datasets considered in the [supplementary materials](https://journals.aps.org/pre/supplemental/10.1103/PhysRevE.104.064136) of *Ellerhoff and Rehfeld (2021)*. 
We thank the research groups for producing and making available their data from model outputs, measurements, paleoclimate, and forcing reconstructions.

Model simulations were obtained from:

- **Y. Zhong et al.**, *Asymmetric Cooling of the Atlantic and Pacific Arctic During the Past Two Millennia: A Dual Observation-Modeling Study*, Geophysical Research Letters (2018)
- **B. L. Otto-Bliesner et al.**, *Climate variability and change since 850 ce an ensemble approach with the community earth system model*, Bulletin of the American Meteorological Society (2016)
-  **J. H. Jungclaus et al.**, *Climate and carbon-cycle variability over the last millennium*, Climate of the Past (2010)
- **J. C. Bühler et al.**, *Comparison of the oxygen isotope signatures in speleothem records and iHadCM3 model simulations for the last millennium*, Climate of the Past (2021)
- **P. Braconnot et al.**, *Strengths and challenges for transient Mid76 to Late Holocene simulations with dynamical vegetation*, Climate of the Past (2019)
- **Z. Liu**, *Transient simulation of last deglaciation with a new mechanism for Bølling–Allerød warming*, Science (2009)
- **N. Fischer and J. H. Jungclaus**, *Evolution of the seasonal temperature cycle in a transient Holocene simulation: Orbital forcing and sea-ice*, Climate of the Past (2011)

Paleoclimate, observational and reanalysis data:

- **H. Hersbach et al.**, *The ERA5 global reanalysis*, Quarterly Journal of the Royal Meteorological Society (2020)
- **C. P. Morice et al.**, *Quantifying uncertainties in global and regional temperature change using an ensemble of observational estimates: The HadCRUT4 data set*, Journal of Geophysical Research Atmospheres (2012)
- **PAGES 2k Consortium**, *Consistent multidecadal variability in global temperature reconstruc50 tions and simulations over the Common Era*, Nature Geoscience (2019)

Forcing reconstructions:
-  **G. A. Schmidt et al.**, *Climate forcing reconstructions for use in PMIP simulations of the Last Millennium (v1.1)*, Geoscientific Model Development (2012)
-  **T. J. Crowley and M. B. Unterman**, *Technical details concerning development of a 1200 yr proxy index for global volcanism*, Earth System Science Data (2013)
- **C. Gao et al.**, *Volcanic forcing of climate over the past 1500 years: An improved ice core-based index for climate models*, Journal of Geophysical Research: At98 mospheres 113, 10.1029/2008JD010239 (2008).
- **M. Toohey and M. Sigl**, *Reconstructed volcanic stratospheric sulfur injections and aerosol optical depth, 500 BCE to 1900 CE, version 2*, World Data Center for Climate (WDCC) at DKRZ (2017)
- **G. Delaygue and E. Bard**, *An Antarctic view of Beryllium-10 and solar activity for the past millennium*, Climate Dynamics (2011)
- **R. Muscheler et al.**, *Solar activity during the last 1000yr inferred from radionuclide records*, Quaternary Science Reviews (2007)
- **F. Steinhilber et al.**, *Total solar irradiance during the Holocene*, Geophysical Research Letters (2009)
- **L. E. Vieira and S. K. Solanki**, *Evolution of the solar magnetic flux on time scales of years to millenia*, Astronomy and Astrophysics (2010)
- **N. A. Krivova et al.**, *Reconstruction of solar total irradiance since 1700 from the surface magnetic flux*, Astronomy and Astrophysics (2007)
- **Y. Wang et al.**, *Modeling the Sun’s Magnetic Field and Irradiance since 1713*, The Astrophysical Journal (2005)
- **C. Fröhlich**, *Solar irradiance variability since 1978: Revision of the PMOD composite during solar cycle 21*, Space Science Reviews (2006)
- **C. D. Keeling et al.**, *Atmospheric carbon dioxide variations at Mauna Loa Observatory*, Hawaii, TELLUS (1976)
- **A. L. Berger**, *Long-term variations of daily insolation and Quaternary climatic changes*, Journal of Atmospheric Sciences (1978)
- To numerically compute orbital variations based on *Berger (1978)*, we used the `Palinsol` package from M. Crucifix (2016)

---

We acknowledge the [R Core team](https://www.R-project.org/) and all package developers of packages used in this study. We thank them for their time and dedication to provide R and the packages to the public. Please see `citation()` for details on the R Core Team and `citation("pkgname")` for details on the developers of individual packages.

The study *Ellerhoff and Rehfed (2021)* has been funded by the [Heidelberg Graduate School for Physics](https://hgsfp.uni-heidelberg.de/), by the [PalMod](https://www.palmod.de/) project (subProject no. 01LP1926C), and by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), [Project No. 395588486](https://gepris.dfg.de/gepris/projekt/395588486?context=projekt&task=showDetail&id=395588486&). It benefited from discussions within the [CVAS](https://pastglobalchanges.org/science/wg/cvas/intro) working group, a working group of the [Past Global Changes (PAGES)](https://pastglobalchanges.org/pal) project. We also thank the members of the [Earth's climate and environmental dynamics](https://www.iup.uni-heidelberg.de/en/research/paleoclimate-dynamics) and the [Climatology and the Biosphere](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/geowissenschaften/arbeitsgruppen/mineralogie-geodynamik/forschungsbereich/klimatologie-und-biosphaere/arbeitsgruppe/) working groups for their helpful tips on preparing the analysis, code, and manuscript. 

Please report bugs to beatrice.ellerhoff(at)iup.uni-heidelberg.de.

*The authors, December 2021*