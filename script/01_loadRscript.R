#written by Dan, Deus and Chikondi
#===============================================================================

#load a package "pacman" used for for installing and loading other packages
if(!require(pacman)) install.packages("pacman")

#load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "reshape2", "dplyr", "tidyr", "minqa", "Metrics", "patchwork", "here",
                        "lmerTest", "lme4", "ggpubr", "janitor", "msm", "rio", "data.table", "png", "grid", "PropCIs"))


#reproduce entire session
addTaskCallback(function(...) {set.seed(1988); TRUE})

#turn off the task call to reset seed if needed
#removeTaskCallback(1)

#fit serotype specific model
source(here("script", "02_st_modelfit.R"))

#fit total carriage model
source(here("script", "03_carr_modelfit.R"))

#describe the study population
source(here("script", "04_carr_descript.R"))

#compare virus circulation in YNHH and study
source(here("script", "05_viral_descript.R"))
