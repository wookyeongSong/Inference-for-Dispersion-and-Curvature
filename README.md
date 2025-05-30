# Inference-for-Dispersion-and-Curvature


## General Information

This repository contains the implementation for the paper **"Inference for Dispersion and Curvature of Random Objects"** which were run using R 4.3.1. It provides two one-click *wrapper* R scripts that reproduce every table and figure in the main manuscript and the supplement.

## Wrapper R scripts

It contains two wrapper R scripts:
*	`DC_wrapper_main.R`: Executes all analyses presented in the main paper and reproduces all the Tables 1-2 and corresponding Figures 2-7, which are saved as PDF files.
*	`DC_wrapper_supplement.R`: Executes all additional simulations and real data analyses in the Supplement of the paper and reproduces the Table S.1 and the corresponding Figures S.1-S.10, which are saved as PDF files.

For running the wrapper R scripts, 
1. Clone or download the repository and set your R working directory to the project root.
2. Download the gait synchronization dataset from [link][https://github.com/deepcharles/gait-data] as referenced in Truong et al. (2019), "A data set for the study of human locomotion with inertial measurement units," Image Processing On Line, 9, 381–390.
3. Set the file path to the downloaded gait data in `Applications/DC_gait.R` file before running the script.
4. Run the wrapper scripts using either `source("DC_wrapper_main.R")` in an R console or `Rscript DC_wrapper_main.R` in a terminal.

## Details on reproducing results, figures, and tables 

In `simulations` folder:
* `/DC_sphere_contam.R`: Produces results corresponding to Figure 2 in Section 4.2.
* `/DC_distribution.R`: Produces results corresponding to Figure 3, Figure 4, Figure S.6, and Figure S.8 in Section 5.1, Section S.5, and Section S.6.
* `/DC_pointcloud.R`: Produces results corresponding to Table 1 and Figure 5 in Section 5.2.
* `/DC_spd.R`: Produces results corresponding to Figure S.1 in Section S.1.1.
* `/DC_sphere.R`: Produces results corresponding to Figure S.2 in Section S.1.2.
* `/DC_power_analysis.R`: Produces results corresponding to Figure S.3 in Section S.1.3.
* `/DC_highdim.R`: Produces results corresponding to Figure S.4 in Section S.1.4.

In `Applications` folder:
* `/DC_gait.R`: Before running this script, set the file path to the gait synchronization data, available at [https://github.com/deepcharles/gait-data]. Running this R code will produce results corresponding to Table 2, Figure 6, Figure S.7 in Section 6.1 and Section S.7.
* `/DC_energy.R`: Produces results corresponding to Figure 7 in Section 6.2.
* `/DC_mnist.R`: Produces results corresponding to Table S.1 in Section S.2.1.
* `/DC_temperature.R`: Produces results corresponding to Figure S.5 in Section S.2.2.

In `src` folder: 
* `/DC_mainfunctions.R`: Contains all main functions necessary for the simulations and real data analysis.

In `Figures` folder: Contains all figures from the main manuscript and supplement saved as `.png` or `.pdf` files. All figures except Figures 1 depend on the `Figures_code` folder.

In `Figures_code` folder: This folder contains R scripts used to reproduce each figure in both the main manuscript and supplement.
* `/Figure2.R`-`/Figure7.R`: Generate Figure 2 through Figure 7 in main manuscript. The figures are automatically saved as PDF files.
* `/FigureS1.R`-`/FigureS10.R`: Generate Figure S.1 through Figure S.10 in the Supplement. The figures are automatically saved as PDF files.
