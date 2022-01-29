IMPORTANT NOTE:

The raw data will be soon available on [Zenodo](https://www.zenodo.org/), together with the code to generate the pre-processed data. 

We kindly ask the user of this repository to download the following datasets into the corresponding directories:

- HadCRUT4: Please download the "HadCRUT.4.6.0.0.anomalies.1.nc" dataset from https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html into the `./processing/raw_data` directory. Our code will call this file as follows: `processing/raw_data/HadCRUT.4.6.0.0.anomalies.1.nc`. The corresponding reference is: - **C. P. Morice et al.**, *Quantifying uncertainties in global and regional temperature change using an ensemble of observational estimates: The HadCRUT4 data set*, Journal of Geophysical Research Atmospheres (2012)
- PAGES2k: 
    - Please download the "PAGES2k_v2.0.0.Rds" from https://springernature.figshare.com/collections/A_global_multiproxy_database_for_temperature_reconstructions_of_the_Common_Era/3285353 into the same directory (`./processing/raw_data`). Our code calls this dataset as `./processing/raw_data/PAGES2k_v2.0.0.Rds`.
    - Please download the reconstruction files ("BHM.txt", "CPS.txt", "DA.txt", "M08.txt", "OIE.txt", "PAI.txt", "PCR.txt") from the "Reconstructions ensemble" folder at https://figshare.com/collections/Global_mean_temperature_reconstructions_over_the_Common_Era/4507043 into the `./processing/raw_data/pages2k/` directory.
    The corresponding references are: **PAGES 2k Consortium**, *Consistent multidecadal variability in global temperature reconstruc50 tions and simulations over the Common Era*, Nature Geoscience (2019) and **PAGES 2k Consortium**, * A global multiproxy database for temperature reconstructions of the Common Era*, Scientific Data (2017)