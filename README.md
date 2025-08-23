# A Statistical Model to Evaluate Radioisotope Contamination Data

This page offers an R program tailored for survival analysis of radioisotope data. The method features:

- rigorous handling of detection limits via cumulative probability distributions,

- extension of survival models to two-component mixtures,

- incorporation of correlations between cesium-134 and cesium-137 using random-effects structures, and

- separation of physical and biological processes in ecological half-life estimation.

We also provide a generalizable formula for quantitative risk evaluation under detection limits.

Note (Windows users):
The TMB package requires compilation of C++ code. On Windows, you need to install Rtools before using this program.

- Download: CRAN Rtools page

- After installation, restart R and ensure that Rtools is on your system PATH (https://cran.r-project.org/bin/windows/Rtools/rtools40.html). <br> Without Rtools, the installation or execution of TMB models may fail.
