# A novel approach to evaluate the transitivity assumption 

## Description of the repository

The repository offers the typical structure of separate folders for data, and R (code/scripts).
* The _data_ folder includes two input text files: 15106232_Taylor(2009).txt and 19637942_Baker(2009).txt;
* The _models_ folder includes two text files, one for the Bayesian random-effects meta-analysis model of continuous outcome and one for the Bayesian random-effects network meta-analysis model of a binary outcome;
* The _R_ folder includes two analysis scripts (Bayesian random-effects MA_Reproducible.R and Bayesian random-effects NMA_Reproducible.R), which source the model and data scripts and perform all analyses, and the five scripts with functions to produce the relevant output; namely, the enhanced balloon-plots, the heatmap of robustness, and the bar-plots with the Kullback-Leibler divergence measure for each scenario analysis.

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the functions of our novel decision framework for robustness of the primary analysis results with emphasis on the summary treatment effects.

## Output 
Prerequisite R packages: [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggpubr](https://cran.r-project.org/web/packages/ggpubr/) and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
