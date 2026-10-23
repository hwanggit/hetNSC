# Simulation studies with correlated informative genes (Sahtout and Wang),
# rewritten with the hetNSC package:
#
#   Example 4 (homoscedastic, both classes correlated, a = 9..23)
#     = inst/legacy/simu_homo_correlated.R
#   Example 5 (heteroscedastic, only class B correlated, a = 9..20)
#     = inst/legacy/simu_heter_correlated.R
#
# Each `iter` is one simulated dataset (the paper uses 100); each dataset is
# analysed with 100 runs per classifier. This is computationally heavy -- run
# different `iter` values as separate jobs:
#
#   Rscript simulation_correlated.R 5 1     # Example 5, iter = 1
#   Rscript simulation_correlated.R 4 1     # Example 4, iter = 1

library(hetNSC)

args    <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) as.integer(args[1]) else 5
iter    <- if (length(args) >= 2) as.integer(args[2]) else 1
stopifnot(example %in% c(4, 5))

N_RUNS     <- 100
CORRELATED <- if (example == 4) "both" else "B"
As         <- if (example == 4) 9:23 else 9:20
OUT_DIR    <- sprintf("hetNSC_example%d_results", example)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

p <- 10000
s <- rep_len(1:20, p)          # gene standard deviations 1, 2, ..., 20 repeated

res <- vector("list", length(As))
for (j in seq_along(As)) {
  a <- As[j]
  set.seed(iter * 8900)
  dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 50, nte2 = 50,
                              ngene = 20, p = p, s = s, a = a,
                              correlated = CORRELATED)
  res[[j]] <- compare_nsc(dat$xtr, dat$ytr, dat$xte, dat$yte,
                          n.runs = N_RUNS,
                          file = file.path(OUT_DIR, sprintf(
                            "simres_a%d_iter%d.rds", a, iter)))
}

# mean test error per method and signal strength
error_tab <- do.call(rbind, lapply(res, function(r) as.data.frame(r$error)))
error_tab$a <- As
error_tab$iter <- iter

# balanced accuracy (%) per method, from the run-averaged confusion matrix
bal_tab <- do.call(rbind, lapply(res, function(r) {
  sm <- summary(r)
  as.data.frame(as.list(stats::setNames(round(100 * sm$balanced.accuracy, 2),
                                        paste0("BalAcc_", sm$method))))
}))
bal_tab$a <- As
bal_tab$iter <- iter

# per-class recall
recall_tab <- do.call(rbind, lapply(seq_along(res), function(j)
  cbind(recall_table(res[[j]]), a = As[j], iter = iter)))

write.csv(error_tab, file.path(OUT_DIR, sprintf("iter%d_error.csv", iter)),
          row.names = FALSE)
write.csv(bal_tab, file.path(OUT_DIR, sprintf("iter%d_balanced_accuracy.csv", iter)),
          row.names = FALSE)
write.csv(recall_tab, file.path(OUT_DIR, sprintf("iter%d_class_sensitivity.csv", iter)),
          row.names = FALSE)
print(error_tab)
print(bal_tab)
