#' Confusion matrix, per-class recall and balanced accuracy
#'
#' `make_cm()` builds a K x K confusion matrix (rows = true class, columns =
#' predicted class). `per_class_recall()` returns the diagonal divided by the
#' row totals (sensitivity of each class) and `balanced_acc()` their mean
#' (macro-averaged recall).
#'
#' @param true_labels,pred_labels Vectors of true and predicted labels.
#' @param classes Class labels defining the rows/columns. Labels not in
#'   `classes` are ignored.
#' @param cm A confusion matrix as returned by `make_cm()`.
#' @return `make_cm()`: an integer matrix; `per_class_recall()`: a named
#'   numeric vector; `balanced_acc()`: a number.
#' @examples
#' cm <- make_cm(c("A", "A", "B", "B"), c("A", "B", "B", "B"), c("A", "B"))
#' cm
#' per_class_recall(cm)
#' balanced_acc(cm)
#' @export
make_cm <- function(true_labels, pred_labels,
                    classes = sort(unique(as.character(true_labels)))) {
  t <- as.character(true_labels)
  p <- as.character(pred_labels)
  ok <- t %in% classes & p %in% classes
  cm <- table(true = factor(t[ok], levels = classes),
              pred = factor(p[ok], levels = classes))
  cm <- unclass(cm)
  storage.mode(cm) <- "integer"
  cm
}

#' @rdname make_cm
#' @export
per_class_recall <- function(cm) diag(cm) / pmax(rowSums(cm), 1L)

#' @rdname make_cm
#' @export
balanced_acc <- function(cm) mean(per_class_recall(cm))

#' Empirical relative heterogeneity rho_k^2 and the crossover condition
#'
#' Computes, for each class k, the empirical relative heterogeneity
#' rho_k^2 = (sigma_k^2 / sigma_bar^2 - 1)^2 of Section 2.4 of the paper,
#' where sigma_bar^2 = sum_k n_k sigma_k^2 / n, summarised across genes, and
#' compares it with the simple crossover bound 2 / (n_k - 1) (obtained with
#' w_k = 1 in eq. 2.11). When `crossover_met` is `TRUE` the class-specific
#' variance estimate is expected to have lower mean squared error than the
#' pooled one, i.e. the heteroscedastic classifiers are expected to help.
#'
#' @param X_train Genes x samples matrix of training data.
#' @param y_train Training class labels (length `ncol(X_train)`).
#' @param X_test,y_test Optional separate data used only to estimate the
#'   class variances (useful when training classes are very small). The
#'   weights n_k always come from the training labels.
#' @param summary_fn Function used to summarise across genes (default
#'   [stats::median()]).
#' @return A data frame with one row per class: `class`, `n_k`, `sigma2_hat`,
#'   `sigma_bar2`, `rho2_hat`, `threshold` (= 2/(n_k - 1)) and
#'   `crossover_met`.
#' @examples
#' set.seed(1)
#' dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 5, nte2 = 5,
#'                             p = 500, a = 12)
#' compute_rho2(dat$xtr, dat$ytr)
#' @export
compute_rho2 <- function(X_train, y_train, X_test = X_train,
                         y_test = y_train, summary_fn = stats::median) {
  stopifnot(is.matrix(X_train), is.matrix(X_test))
  stopifnot(length(y_train) == ncol(X_train))
  stopifnot(length(y_test) == ncol(X_test))
  stopifnot(nrow(X_train) == nrow(X_test))

  y_train <- as.character(y_train)
  y_test <- as.character(y_test)
  classes <- sort(unique(y_train))
  p <- nrow(X_train)

  n_k <- vapply(classes, function(k) sum(y_train == k), numeric(1))
  n <- sum(n_k)

  s2 <- vapply(classes, function(k) {
    Xk <- X_test[, y_test == k, drop = FALSE]
    if (ncol(Xk) < 2L) {
      warning(sprintf(paste0("Class '%s' has only %d sample(s) in the ",
                             "variance-estimation data; its variance ",
                             "estimates are NA."), k, ncol(Xk)),
              call. = FALSE)
      return(rep(NA_real_, p))
    }
    apply(Xk, 1L, stats::var)
  }, numeric(p))
  s2 <- matrix(s2, nrow = p, dimnames = list(NULL, classes))

  sigma_bar2_vec <- as.numeric(s2 %*% n_k) / n
  rho2_mat <- (sweep(s2, 1L, sigma_bar2_vec, "/") - 1)^2

  rho2_summary <- apply(rho2_mat, 2L, summary_fn, na.rm = TRUE)
  sigma2_summary <- apply(s2, 2L, summary_fn, na.rm = TRUE)
  sigma_bar2_val <- summary_fn(sigma_bar2_vec, na.rm = TRUE)
  threshold <- ifelse(n_k > 1, 2 / (n_k - 1), Inf)

  data.frame(
    class = classes,
    n_k = as.integer(n_k),
    sigma2_hat = round(unname(sigma2_summary), 4),
    sigma_bar2 = round(sigma_bar2_val, 4),
    rho2_hat = round(unname(rho2_summary), 4),
    threshold = round(unname(threshold), 4),
    crossover_met = unname(rho2_summary > threshold),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
