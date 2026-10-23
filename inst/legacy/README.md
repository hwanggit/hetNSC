# Original research code

These are the original scripts used for the paper, kept unchanged for
reference and reproducibility. They are **not** loaded by the package.

| File | Contents |
|------|----------|
| pamr.soft.R, pamr.hard.R, pamr.order.R | Homoscedastic NSC with soft / hard / order thresholding and deep CV (STh, HTh, OTh) |
| pamr.soft.hetero.R, pamr.hard.hetero.R, pamr.order.hetero.R | Heteroscedastic NSC (STh3, HTh3, OTh3) |
| simu_heter_correlated.R, sahtout1_run_simu_corr.R | Correlated-gene simulation, class B correlated (heteroscedastic) |
| simu_homo_correlated.R | Correlated-gene simulation, both classes correlated (homoscedastic) |

Mapping to the package:

| Original | hetNSC |
|----------|--------|
| `my.pamr.Scv.old / .hard / .order / .heter / .H.heter / .O.heter` | `hnsc_select(x, y, method = "STh" / "HTh" / "OTh" / "STh3" / "HTh3" / "OTh3")` |
| `pamr.train.*` | `hnsc_train()` |
| `pamr.cv.*` | `hnsc_cv()` |
| `pamr.predict*` | `predict()` on an `hnsc` object |
| `generate_sample()` | `simulate_hetero_data(correlated = "B")` or `"both"` (returns genes x samples) |
| `run_one_dataset()` | `compare_nsc()` |
| `mean_cm()` + `build_recall_table()` / `build_balanced_acc_table()` | `recall_table()`, `summary()` of a `compare_nsc()` result |
| `compute_rho2()`, `make_cm()`, `per_class_recall()`, `balanced_acc()` | same names |

The package test `tests/testthat/test-legacy.R` checks that the package
reproduces these scripts exactly (same thresholds, CV errors and posterior
probabilities for the same random seed).
