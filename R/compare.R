#' Compare homoscedastic and heteroscedastic NSC classifiers on one dataset
#'
#' Reproduces the evaluation protocol of the paper: for each thresholding
#' rule, the homoscedastic classifier (STh/HTh/OTh) and its heteroscedastic
#' counterpart (STh3/HTh3/OTh3) are fitted `n.runs` times, each time with the
#' threshold chosen by [hnsc_select()] after a random reordering of the
#' training samples (which changes the cross-validation folds). Each fit is
#' evaluated on the test set.
#'
#' @param x,y Training data (genes x samples) and labels.
#' @param xtest,ytest Test data (genes x samples) and labels.
#' @param n.runs Number of repetitions (100 in the paper).
#' @param thresholding Which thresholding rules to run, in this order.
#' @param seed_fun Function of the run index returning the seed set before
#'   that run. The default, `function(i) 1000 + 10 * i`, is the one used in
#'   the paper. Use `NULL` to not set seeds.
#' @param file Optional path; if given the result is also saved with
#'   [saveRDS()].
#' @param verbose Print progress.
#' @return A list of class `"hnsc_comparison"` with
#'   * `error`: mean test error rate of each method,
#'   * `ave_num_genes`: mean number of genes used,
#'   * `errors`, `n_genes`: the per-run values (runs x methods matrices),
#'   * `cm_list`: per method, a list of per-run confusion matrices,
#'   * `classes`, `n` (training sample size) and `K`.
#'
#'   Use [summary()] on the result for a table that includes balanced
#'   accuracy, and [recall_table()] for per-class recall.
#' @details The training samples are re-shuffled cumulatively across runs and
#'   thresholding rules, exactly as in the script used for the paper, so that
#'   `compare_nsc()` reproduces those numbers.
#' @examples
#' \donttest{
#' set.seed(8900)
#' dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 50, nte2 = 50,
#'                             p = 300, a = 15)
#' res <- compare_nsc(dat$xtr, dat$ytr, dat$xte, dat$yte, n.runs = 2,
#'                    verbose = FALSE)
#' summary(res)
#' }
#' @export
compare_nsc <- function(x, y, xtest, ytest, n.runs = 100,
                        thresholding = c("hard", "soft", "order"),
                        seed_fun = function(i) 1000 + 10 * i,
                        file = NULL, verbose = TRUE) {
  x <- check_xy(x, y)
  if (is.data.frame(xtest)) xtest <- as.matrix(xtest)
  if (nrow(xtest) != nrow(x)) stop("`xtest` must have the same genes (rows) as `x`.",
                                   call. = FALSE)
  if (length(ytest) != ncol(xtest)) stop("length(ytest) must equal ncol(xtest).",
                                        call. = FALSE)
  thresholding <- match.arg(thresholding, several.ok = TRUE)
  cl <- sort(unique(as.character(y)))
  ytest <- as.character(ytest)
  labels <- as.vector(rbind(
    vapply(thresholding, method_label, "", variance = "homo"),
    vapply(thresholding, method_label, "", variance = "hetero")))
  errors <- n_genes <- matrix(NA_real_, n.runs, length(labels),
                              dimnames = list(NULL, labels))
  cm_list <- stats::setNames(lapply(labels, function(m) vector("list", n.runs)),
                             labels)
  y <- as.character(y)

  for (meth in thresholding) {
    for (itr in seq_len(n.runs)) {
      if (!is.null(seed_fun)) set.seed(seed_fun(itr))
      ord <- sample(length(y))
      x <- x[, ord, drop = FALSE]
      y <- y[ord]
      yf <- factor(y)
      for (v in c("homo", "hetero")) {
        lab <- method_label(meth, v)
        sel <- deep_cv(x, yf, meth, v)
        fit <- hnsc_train(x, y, method = meth, variance = v,
                          threshold = sel$threshold)
        pred <- as.character(predict.hnsc(fit, xtest, threshold = sel$threshold))
        cm_list[[lab]][[itr]] <- make_cm(ytest, pred, cl)
        n_genes[itr, lab] <- fit$nonzero
        errors[itr, lab] <- mean(pred != ytest)
      }
      if (verbose)
        message(sprintf("%-5s run %d/%d: %s", meth, itr, n.runs,
                        paste(sprintf("%s=%.3f", colnames(errors),
                                      errors[itr, ])[!is.na(errors[itr, ])],
                              collapse = " ")))
    }
  }

  res <- list(error = as.list(colMeans(errors)),
              ave_num_genes = as.list(colMeans(n_genes)),
              errors = errors, n_genes = n_genes, cm_list = cm_list,
              classes = cl, n = ncol(x), K = length(cl))
  class(res) <- "hnsc_comparison"
  if (!is.null(file)) saveRDS(res, file)
  res
}

# Average confusion matrix over runs (as in the paper's summary code).
mean_cm <- function(cm_runs) {
  valid <- Filter(Negate(is.null), cm_runs)
  if (!length(valid)) return(NULL)
  Reduce("+", valid) / length(valid)
}

#' Summaries of a classifier comparison
#'
#' `summary()` gives, for each method, the mean and standard deviation of the
#' test error over runs, the balanced accuracy and the mean number of genes
#' used. `recall_table()` gives the recall (sensitivity) of every class.
#' Both are computed from the confusion matrix averaged over runs, as in the
#' paper.
#'
#' @param object A result of [compare_nsc()].
#' @param ... Not used.
#' @return A data frame. `summary()`: one row per method. `recall_table()`:
#'   one row per class and one `recall_<method>` column per method.
#' @examples
#' \donttest{
#' set.seed(8900)
#' dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 50, nte2 = 50,
#'                             p = 300, a = 15)
#' res <- compare_nsc(dat$xtr, dat$ytr, dat$xte, dat$yte, n.runs = 2,
#'                    thresholding = "hard", verbose = FALSE)
#' summary(res)
#' recall_table(res)
#' }
#' @export
summary.hnsc_comparison <- function(object, ...) {
  labs <- colnames(object$errors)
  bal <- vapply(labs, function(m) balanced_acc(mean_cm(object$cm_list[[m]])),
                numeric(1))
  data.frame(method = labs,
             error = colMeans(object$errors),
             error.sd = apply(object$errors, 2, stats::sd),
             balanced.accuracy = unname(bal),
             genes = colMeans(object$n_genes),
             row.names = NULL)
}

#' @rdname summary.hnsc_comparison
#' @export
recall_table <- function(object) {
  stopifnot(inherits(object, "hnsc_comparison"))
  out <- data.frame(class = object$classes, stringsAsFactors = FALSE)
  for (m in names(object$cm_list))
    out[[paste0("recall_", m)]] <-
      unname(per_class_recall(mean_cm(object$cm_list[[m]]))[object$classes])
  out
}

#' @export
print.hnsc_comparison <- function(x, ...) {
  cat(sprintf("NSC comparison: %d runs, %d classes, %d training samples\n",
              nrow(x$errors), x$K, x$n))
  print(summary(x), row.names = FALSE, digits = 4)
  invisible(x)
}
