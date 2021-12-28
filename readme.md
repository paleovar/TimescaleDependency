Coming soon. 

This repository provides code to reproduce the figures of "Probing the timescale dependency of local and global variations in surface air temperature from climate simulations and reconstructions of the last millennia" accepted in Physical Review E (2021) (preprint https://arxiv.org/abs/2102.00458).

Authors: Beatrice Ellerhoff and Kira Rehfeld 

Contact: beatrice.ellerhoff(at)iup.uni-heidelberg.de


required packages:
- dplyr
- tibble
- tidyr
- zoo
- RColorBrewer
- ggplot2
- latex2exp
- tseries
- PaleoSpec

Creating some of the figures requires:
- cowplot
- forcats
- ggnewscale
- rgdal (in case you encounter difficulties in the installation, this post might help: https://gist.github.com/dncgst/111b74066eaea87c92cdc5211949cd1e)
- sp
- sf (in case you encounter difficulties, this link might be useful: https://r-spatial.github.io/sf/)
- viridisLite
- raster
- ggcorrplot
- nest