# BARDSEM
## Zachary Roman <ZacharyJoseph.Roman@uzh.ch>

This repo houses the model files and run script from the homicide example in Roman and Brandt (2020) for technical details please see the paper at https://doi.org/10.1080/00273171.2021.1957663

These files pre-proccess and estimate a Bayesian Auto-regressive Dependence Structural Equation Model (BARDSEM).

PreProcRunFile.R - This is what you should open to process the data and run the model,
there are also cloro-plots from the paper in this script.

SEM_ENDLAG_impacts.stan - This is the model syntax.

data/shapeData - This folder contains the shape files originally developed for GIS software
R imports these fine, I use these to develop the spatial relationship matricies in the paper.

UCR_arrests - This is the Uniform Crime Reporting (UCR) data which I use to obtain regional crime rates
The codebook is in this folder which breaks down reporting standards and variable definitions. 

homocide1990.dta12 - This is the Homocide data file, I end up not using it int he paper, but there
are some interesting plots with this data that I considered adding to the paper. I eventualy
chose not to include this. 
