# -------------------------------------------------------------------------
# Title : Supplement of ``Inference for Dispersion and Curvature of Random Objects"
# Description : Wrapper file that sequentially runs through each section of the analyses done and produce the figures in the Supplement
# -------------------------------------------------------------------------

## Import main functions and packages
source("src/DC_mainfunctions.R")


# -------------------------------------------------------------------------
# Additional simulation results: Spaces of symmetric positive definite matrices
#
# Manuscript reference : Section S.1.1
# Figure reproduced    : Figure S.1
# Simulation code reference : simulations/DC_spd.R
# -------------------------------------------------------------------------

## Conduct ambient curvature test for spaces of symmetric positive definite matrices with 6 different distances in Section S.1.1
source("simulations/DC_spd.R") 
source("Figures_code/FigureS1.R") # Generate Figure S.1


# -------------------------------------------------------------------------
# Additional simulation results: Spherical data with geodesic distance
#
# Manuscript reference : Section S.1.2
# Figure reproduced    : Figure S.2
# Simulation code reference : simulations/DC_sphere.R
# -------------------------------------------------------------------------

## Conduct ambient curvature test for Spherical data with geodesic distance in Section S.1.2
source("simulations/DC_sphere.R") 
source("Figures_code/FigureS2.R") # Generate Figure S.2


# -------------------------------------------------------------------------
# Power and Type-I error analysis under model spaces
#
# Manuscript reference : Section S.1.3
# Figure reproduced    : Figure S.3
# Simulation code reference : simulations/DC_power_analysis.R
# -------------------------------------------------------------------------

## Conduct power and Type-I error analysis under model spaces in Section S.1.3
source("simulations/DC_power_analysis.R") 
source("Figures_code/FigureS3.R") # Generate Figure S.3


# -------------------------------------------------------------------------
# Power and Type-I error analysis under High-dimensional ambient Euclidean random objects with low intrinsic dimension
#
# Manuscript reference : Section S.1.4
# Figure reproduced    : Figure S.4
# Simulation code reference : simulations/DC_highdim.R
# -------------------------------------------------------------------------

## Conduct power and Type-I error analysis under High-dimensional ambient Euclidean random objects with low intrinsic dimension in Section S.1.4
source("simulations/DC_highdim.R") 
source("Figures_code/FigureS4.R") # Generate Figure S.4


# -------------------------------------------------------------------------
# Additional real data analysis: MNIST Dataset
#
# Data available at R package `dslabs`
# Manuscript reference : Section S.2.1
# Table reproduced     : Table S.1
# Data application code reference : Applications/DC_mnist.R
# -------------------------------------------------------------------------

## MNIST data analysis in Section S.2.1
source("Applications/DC_mnist.R") # Print Table S.1


# -------------------------------------------------------------------------
# Additional real data analysis: US airport weather station temperature data
#
# Data available at Applications/airport data
# Manuscript reference : Section S.2.2
# Figure reproduced    : Figure S.5
# Data application code reference : Applications/DC_temperature.R
# -------------------------------------------------------------------------

## US airport weather station temperature data analysis in Section S.2.2
source("Applications/DC_temperature.R")
source("Figures_code/FigureS5.R") # Generate Figure S.5


# -------------------------------------------------------------------------
# Verification of conditions for distributional simulations
#
# Manuscript reference : Section S.5
# Figure reproduced    : Figure S.6
# Simulation code reference : simulations/DC_distribution.R
# -----------------------------------------------------------------------

## ISOMAP representation plot in Section S.5
source("simulations/DC_distribution.R") 
source("Figures_code/FigureS6.R") # Generate Figure S.6


# -------------------------------------------------------------------------
# Sensitivity analysis for the choice of input distance
#
# Data available at https://github.com/deepcharles/gait-data
# Reference: Truong et al. 2019, 'A data set for the study of human locomotion with inertial measurements units', Image Processing On Line 9, 381–390."
# Manuscript reference : Section S.6
# Figures reproduced    : Figure S.7, Figure S.8
# Codes reference : Applications/DC_gait.R, simulations/DC_distribution.R
# -------------------------------------------------------------------------

## Sensitivity analysis for the choice of input distance using Gait synchronization data in the space of symmetric positive definite matrices
source("Applications/DC_gait.R")
source("Figures_code/FigureS7.R") # Generate Figure S.7

## Sensitivity analysis for the choice of input distance using distributional simulations with Wasserstein distance
source("simulations/DC_distribution.R") 
source("Figures_code/FigureS8.R") # Generate Figure S.8


# -------------------------------------------------------------------------
# Relation between metric curvature and Alexandrov curvature
#
# Manuscript reference : Section S.7
# Figures reproduced    : Figure S.9, Figure S.10
# -------------------------------------------------------------------------

source("Figures_code/FigureS9.R") # Generate Figure S.9
source("Figures_code/FigureS10.R") # Generate Figure S.10











