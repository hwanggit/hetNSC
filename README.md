# hetNSC
Heteroscedastic and homoscedastic nearest shrunken centroid classifiers with hard, soft, and order thresholding

# R-package-Heteroscedastic-and-homoscedastic-Shrunken-Centroid-Classifiers

# Package name: hetNSC

Heteroscedastic nearest shrunken centroid (NSC) classification for
high-dimensional, low-sample-size data such as gene expression.

The package implements six classifiers:

| Label | Thresholding | Variance in test statistic and discriminant |
|-------|--------------|---------------------------------------------|
| STh   | soft  | pooled: the original PAM (Tibshirani et al. 2002) |
| HTh   | hard  | pooled (Sahtout, Wang and Ghimire 2024) |
| OTh   | order | pooled (Sahtout, Wang and Ghimire 2024) |
| STh3  | soft  | class-specific (Sahtout and Wang) |
| HTh3  | hard  | class-specific (Sahtout and Wang) |
| OTh3  | order | class-specific (Sahtout and Wang) |

The threshold is chosen by the refined ("deep") cross-validation search of
Sahtout, Wang and Ghimire (2024).

## Installation

```r
# from a local copy of the source package
install.packages("hetNSC_0.1.0.tar.gz", repos = NULL, type = "source")

# or, once it is on GitHub
# remotes::install_github("hwanggit/hetNSC")
```

The package only depends on base R and **MASS**.

## Quick start

Expression matrices have **genes in rows and samples in columns** (as in
**pamr**).

```r
library(hetNSC)

set.seed(8900)
dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 50, nte2 = 50,
                            p = 1000, a = 14)

fit <- hnsc(dat$xtr, dat$ytr, method = "HTh3")   # or method = "hard", variance = "hetero"
fit
pred <- predict(fit, dat$xte)                    # class labels
mean(pred != dat$yte)                            # test error
predict(fit, dat$xte, type = "posterior")        # posterior probabilities
fit$selected.genes                               # genes used
```

## Main functions

| Function | Purpose |
|----------|---------|
| `hnsc()` | Fit a classifier with the threshold selected by deep CV |
| `predict()` | Class labels, posterior probabilities, shrunken centroids or selected genes |
| `hnsc_train()`, `hnsc_cv()`, `hnsc_select()` | Threshold grid fit, cross-validation, deep CV search |
| `compare_nsc()` | Repeated-run comparison of homoscedastic vs heteroscedastic versions (protocol of the paper); `summary()` and `recall_table()` summarise it |
| `simulate_hetero_data()` | Correlated-gene simulation designs (class B only, or both classes) |
| `compute_rho2()` | Relative heterogeneity and bias-variance crossover check (Section 2.4) |
| `make_cm()`, `per_class_recall()`, `balanced_acc()` | Class-level metrics |
| `hnsc_methods()` | Table of the six classifiers |

See `vignette("hetNSC")` for a walk-through, and
`system.file("scripts", "simulation_correlated.R", package = "hetNSC")` for the
full correlated-gene simulation studies.

## Reproducibility

The original research scripts are kept in `inst/legacy`. For the same random
seed the package returns the same selected thresholds, cross-validation
errors, test errors, numbers of genes and posterior probabilities as those
scripts; `tests/testthat/test-legacy.R` checks this.

## Citation

```r
citation("hetNSC")
```

* Sahtout, M. O. and Wang, H. Heteroscedastic nearest shrunken centroid
  classification: Methods and applications to high-dimensional data.
  *Statistical Applications in Genetics and Molecular Biology* (accepted).
* Sahtout, M. O., Wang, H. and Ghimire, S. (2024). Different thresholding
  methods on nearest shrunken centroid algorithm. *Communications in
  Statistics - Simulation and Computation*, 53(3), 1444-1460.
  <https://doi.org/10.1080/03610918.2022.2047201>

## License

GPL (>= 2). Parts of the code are adapted from the **pamr** package by
Hastie, Tibshirani, Narasimhan and Chu.

