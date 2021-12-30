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

#graphical parameters for plotting
pointsize=0.4
notationsize=3
#ggplot theme 
theme_td <- function(textsize=9){
  theme_classic(base_size=textsize) +
    theme(axis.title = element_text(size = textsize),
          axis.text = element_text(size = textsize),
          axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"), color="black"),
          axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), color="black"),
          axis.ticks.length = unit(-1.4, "mm")) +
    theme(legend.key = element_rect(color = "transparent"),
          legend.background=element_blank(), 
          legend.title = element_blank(),
          legend.text = element_text(size=textsize),
          legend.key.height = unit(0.2, "cm")) 
}
