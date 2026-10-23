#' hetNSC: Heteroscedastic Nearest Shrunken Centroid Classification
#'
#' Nearest shrunken centroid classifiers for high-dimensional data, with
#' either the pooled-variance (homoscedastic) test statistic and discriminant
#' function of PAM or the heteroscedastic (class-specific variance) ones of
#' Sahtout and Wang, each combined with soft, hard or order thresholding.
#'
#' @section The six classifiers:
#' * **STh** (soft, homoscedastic): the original PAM classifier of
#'   Tibshirani et al. (2002).
#' * **HTh**, **OTh** (hard / order, homoscedastic): the homoscedastic
#'   classifiers of Sahtout, Wang and Ghimire (2024), which replace PAM's
#'   soft thresholding with hard or order thresholding.
#' * **STh3**, **HTh3**, **OTh3**: the heteroscedastic classifiers of
#'   Sahtout and Wang.
#'
#' So `variance = "homo"` is the original PAM only when `method = "soft"`.
#'
#' @section Main functions:
#' * [hnsc()]: fit a classifier, choosing the threshold by refined
#'   cross-validation; then [predict.hnsc()].
#' * [hnsc_train()], [hnsc_cv()], [hnsc_select()]: lower-level steps.
#' * [compare_nsc()]: the repeated-run comparison protocol of the paper;
#'   [summary.hnsc_comparison()] and [recall_table()] summarise it.
#' * [simulate_hetero_data()]: the correlated-gene simulation designs.
#' * [compute_rho2()]: the bias-variance crossover diagnostic of Section 2.4.
#' * [make_cm()], [per_class_recall()], [balanced_acc()]: class-level metrics.
#'
#' @section Data layout:
#' As in \pkg{pamr}, expression matrices have **genes in rows and samples in
#' columns**.
#'
#' @keywords internal
"_PACKAGE"
