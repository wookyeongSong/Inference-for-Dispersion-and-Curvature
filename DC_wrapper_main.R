# -------------------------------------------------------------------------
# Title : ``Inference for Dispersion and Curvature of Random Objects"
# Description : Wrapper file that sequentially runs through each section of the analyses done and produce the figures in the main paper 
# -------------------------------------------------------------------------

## Import main functions and packages
source("src/DC_mainfunctions.R")


# -------------------------------------------------------------------------
# Simulation Study: Noise-Contaminated Spherical Data
#
# Manuscript reference : Section 4.2
# Figure reproduced    : Figure 2
# Simulation code reference : simulations/DC_sphere_contam.R
# -------------------------------------------------------------------------

## Conduct intrinsic curvature test for Noise-Contaminated Spherical Data in Section 4.2
source("simulations/DC_sphere_contam.R") 
source("Figures_code/Figure2.R") # Generate Figure 2


# -------------------------------------------------------------------------
# Simulation Study: Distributional data with 2-Wasserstein distance
#
# Manuscript reference : Section 5.1
# Figures reproduced    : Figure 3, Figure 4
# Simulation code reference : simulations/DC_distribution.R
# -------------------------------------------------------------------------

## Conduct intrinsic curvature test for Distributional data with 2-Wasserstein distance in Section 5.1
source("simulations/DC_distribution.R") 
source("Figures_code/Figure3.R") # Generate Figure 3
source("Figures_code/Figure4.R") # Generate Figure 4


# -------------------------------------------------------------------------
# Simulation Study: Point cloud data with intrinsic distance
#
# Manuscript reference : Section 5.2
# Figure reproduced    : Figure 5
# Table reproduced     : Table 1
# Simulation code reference : simulations/DC_pointcloud.R
# -------------------------------------------------------------------------

## Conduct intrinsic curvature test for Distributional data with Point cloud data with intrinsic distance in Section 5.2
source("simulations/DC_pointcloud.R") # Print Table 1
source("Figures_code/Figure5.R") # Generate Figure 5


# -------------------------------------------------------------------------
# Data Application: Gait synchronization analysis
#
# Data available at https://github.com/deepcharles/gait-data
# Reference: Truong et al. 2019, 'A data set for the study of human locomotion with inertial measurements units', Image Processing On Line 9, 381–390."
# Manuscript reference : Section 6.1
# Figure reproduced    : Figure 6
# Table reproduced     : Table 2
# Data application code reference : Applications/DC_gait.R
# -------------------------------------------------------------------------

## Gait synchronization analysis in Section 6.1
source("Applications/DC_gait.R") # Print Table 2
source("Figures_code/Figure6.R") # Generate Figure 6


# -------------------------------------------------------------------------
# Data Application: Energy source data analysis
#
# Data available at Applications/energy_data
# Manuscript reference : Section 6.2
# Figure reproduced    : Figure 7
# Data application code reference : Applications/DC_energy.R
# -------------------------------------------------------------------------

## Energy source data analysis in Section 6.2
source("Applications/DC_energy.R")
source("Figures_code/Figure7.R") # Generate Figure 7

