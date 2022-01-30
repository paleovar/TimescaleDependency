The raw data will be soon available on [Zenodo](https://www.zenodo.org/), together with the code to generate the pre-processed data. 

**Important note on downloading external datasets prior to data generation:**
We kindly ask the user of this repository to download the following datasets into the corresponding directories:

- HadCRUT4: Please download the "HadCRUT.4.6.0.0.anomalies.1.nc" dataset from https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html into the `./processing/raw_data` directory. Our code will call this file as follows: `processing/raw_data/HadCRUT.4.6.0.0.anomalies.1.nc`. The corresponding reference is: - **C. P. Morice et al.**, *Quantifying uncertainties in global and regional temperature change using an ensemble of observational estimates: The HadCRUT4 data set*, Journal of Geophysical Research Atmospheres (2012)
- PAGES2k: 
    - Please download the "PAGES2k_v2.0.0.Rds" from https://springernature.figshare.com/collections/A_global_multiproxy_database_for_temperature_reconstructions_of_the_Common_Era/3285353 into the same directory (`./processing/raw_data`). Our code calls this dataset as `./processing/raw_data/PAGES2k_v2.0.0.Rds`.
    - Please download the reconstruction files ("BHM.txt", "CPS.txt", "DA.txt", "M08.txt", "OIE.txt", "PAI.txt", "PCR.txt") from the "Reconstructions ensemble" folder at https://figshare.com/collections/Global_mean_temperature_reconstructions_over_the_Common_Era/4507043 into the `./processing/raw_data/pages2k/` directory.
    The corresponding references are: **PAGES 2k Consortium**, *Consistent multidecadal variability in global temperature reconstruc50 tions and simulations over the Common Era*, Nature Geoscience (2019) and **PAGES 2k Consortium**, * A global multiproxy database for temperature reconstructions of the Common Era*, Scientific Data (2017)
    
**Important note on how to generate the data:**
To generate a datase`*.Rds`, run the script `get_*.Rds`. Please be aware that we have equipped the data generating code with a few "safety mechanisms". You will find a `save <- F` at the top of data generating scripts. Please carefully change this to `save <- T` and mind that this will potentially overwrite excisting files. 
Many data generating scripts include "for-loops" over the datasets or iterations (i.e. bootstrapping). Please carefully choose your parameters and try first with single datasets. Some of the codes might require large memory and long computing times. There are warnings in the according scripts. 

**References:**
The data files in `./raw_data/` were obtained by extracting the temperature fields from the following publications. In particular, we alligned the time axis, computed temperature anomalies, cleaned the data (i.e. performed linear interpolation of NAs if necessary) and shaped the data into a consistent format to facilitate the analysis.

- CESM 1 past2k: **Y. Zhong et al.**, *Asymmetric Cooling of the Atlantic and Pacific Arctic During the Past Two Millennia: A Dual Observation-Modeling Study*, Geophysical Research Letters (2018) and 
- CESM-LME: **B. L. Otto-Bliesner et al.**, *Climate variability and change since 850 ce an ensemble approach with the community earth system model*, Bulletin of the American Meteorological Society (2016)
- MPI-ESM ML: **J. H. Jungclaus et al.**, *Climate and carbon-cycle variability over the last millennium*, Climate of the Past (2010)
- HadCM3 LM1: **J. C. Bühler et al.**, *Comparison of the oxygen isotope signatures in speleothem records and iHadCM3 model simulations for the last millennium*, Climate of the Past (2021)
- IPSL-p6k: **P. Braconnot et al.**, *Strengths and challenges for transient Mid76 to Late Holocene simulations with dynamical vegetation*, Climate of the Past (2019)
- TraCE-21k: **Z. Liu**, *Transient simulation of last deglaciation with a new mechanism for Bølling–Allerød warming*, Science (2009)
- ECH5/MPIOM-p6k **N. Fischer and J. H. Jungclaus**, *Evolution of the seasonal temperature cycle in a transient Holocene simulation: Orbital forcing and sea-ice*, Climate of the Past (2011)
- ERA5: **H. Hersbach et al.**, *The ERA5 global reanalysis*, Quarterly Journal of the Royal Meteorological Society (2020)
