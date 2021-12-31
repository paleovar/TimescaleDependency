#load required packages
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(zoo)
library(scales)
library(tseries)
library(PaleoSpec)
library(tibble)
select <- dplyr::select

#load meta data
signal_tbb <- readRDS("helpers/signal_tbb.Rds")

