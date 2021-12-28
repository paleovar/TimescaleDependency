library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(scales)
library(tseries)
library(PaleoSpec)
select <- dplyr::select

signal_tbb <- readRDS("helpers/signal_tbb.Rds")
model <- c(signal_tbb %>% filter(type=="model") %>% select(signal))$signal
obs <- c(signal_tbb %>% filter(type=="obs") %>% select(signal))$signal
signal <- c(model, obs)

#graphical parameters
pointsize=0.4
#textsize=9
notationsize=3

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