# A novel approach to evaluate the transitivity assumption 

## Description of the repository

The repository offers the typical structure of separate folders for data, and R (code/scripts).
* The _data_ folder includes four .RData files: _singh_dataset_, and _singh_network_ for the first motivating example, and _baker_dataset_, and _baker_network_ for the second motivating example
* The _R_ folder includes one function (long.to.wide_function.R) and three analysis R scripts to replicate the Figures (main and supplementary): _A. Main Figures_Singh & Baker_ to replicate the main figures, _B. Supplementary Figures_Singh_ to replicate the supplementary figures of the first example and _C. Supplementary Figures_Baker_ to replicate the supplementary figures of the second example.

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the functions of our novel decision framework for robustness of the primary analysis results with emphasis on the summary treatment effects.

## Output 
Prerequisite R packages: [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggpubr](https://cran.r-project.org/web/packages/ggpubr/) and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
