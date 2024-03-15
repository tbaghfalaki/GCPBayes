# Descripsion 
Execute a Gibbs sampler to analyze a multivariate Bayesian sparse group selection model, employing Dirac, continuous, and hierarchical spike priors to detect pleiotropic effects on two traits. This package is specifically designed for summary statistics comprising estimated regression coefficients and their corresponding covariance matrix.

### Installation
To acquire the latest development version of GCPBayes, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/GCPBayes")
```
This will seamlessly fetch and install the most up-to-date version of GCPBayes for your use.

The package is now accessible on R CRAN as well![https://cran.r-project.org/web/packages/GCPBayes/index.html]
```
  install.packages("GCPBayes")
```

### References 
#### The methodology paper
Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta‐analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. *Statistics in Medicine*, 40(6), 1498-1518. [https://onlinelibrary.wiley.com/doi/10.1002/sim.8855]

#### A Vignette
Baghfalaki, T., Sugier, PE, Asgari, Y., Truong, T., & Liquet, B. (2023). GCPBayes: An R package for studying Cross-Phenotype Genetic Associations with Group-level Bayesian Meta-Analysis. *RJournal* , 15 (1), 122-141. [https://journal.r-project.org/articles/RJ-2023-028/RJ-2023-028.pdf]

#### A Pipeline 
Asgari, Y., Sugier, P. E., Baghfalaki, T., Lucotte, E., Karimi, M., Sedki, M., ... & Truong, T. (2023). GCPBayes pipeline: a tool for exploring pleiotropy at the gene level. *NAR Genomics and Bioinformatics*, 5(3), lqad065. [https://academic.oup.com/nargab/article/5/3/lqad065/7219410]

#### Declaration and Contributions
Taban Baghfalaki has made significant contributions to the development of this package, including the design and implementation of key functionalities such as CS, DS, and HS methods, along with their summaries and associated plots.

As of version 4.2.0 (2024-03-14), Taban Baghfalaki has stepped down as the maintainer of the package.


