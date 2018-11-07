# Simple simulations for comparing scenarios for polygenic risk scores vs IVW MR estimates

Requires

```r
devtools::install_github("explodecomputer/simulateGP")
devtools::install_github("MRCIEU/TwoSampleMR")
install.packages("tidyverse")
```


We ran simulations to explore the differences in performance between the IVW random effects MR model and polygenic risk score analysis under different levels of horizontal pleiotropy. We ran 1000 simulations in which trait X had a causal effect on trait Y, and 1000 null simulations in which there was no causal effect. For each simulation we created a pair of datasets, so that genetic effects could be estimated in the first and then predicted in the second. The parameters for the simulations were as follows.

For each simulation both discovery and test datasets had the same underlying model. Each dataset comprised 10000 samples with data on 50 SNPs, and phenotypes X and Y. Each SNP had a random effect on trait X and all SNPs together explained 15% of the phenotypic variance. The set of 50 SNPs exhibited balanced horizontal pleiotropic effects on trait Y such that

$$
\beta_g_{j}y ~ N(0, \sigma^2_{h})
$$

and for each simulation a random value was sampled for \sigma^2_{h} between 0 and 0.2. The effect SNP effects on X were obtained from the discovery dataset using an additive genetic model. Only SNPs with p < 1e-8 were retained for MR and PRS analysis. For MR analysis the SNP effects on Y were obtained from the test dataset and then two-sample MR estimates were obtained. For PRS analysis, the genetic score was constructed in the test dataset using the SNP-exposure effect estimates from the discovery dataset as weights. All code is available here [https://github.com/explodecomputer/prs-vs-mr](https://github.com/explodecomputer/prs-vs-mr).


## To run

```r
Rscript prs_vs_mr_simulation.run.r
Rscript prs_vs_mr_simulation.plots.r
```

Creates:

```
prs_vs_mr_simulation.rdata
prs_vs_mr_simulation.plot1.pdf
prs_vs_mr_simulation.plot2.pdf
```

