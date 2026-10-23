# Internal helpers. Most are taken (with minimal changes) from the 'pamr'
# package and from the original research code of Sahtout and Wang, so that the
# package reproduces the published results exactly, including the random
# number stream used for cross-validation fold assignment.

# ---- thresholding rules ---------------------------------------------------

soft_shrink <- function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- sign(delta) * dif * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}

hard_shrink <- function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- delta * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}

# Keep the `threshold` largest |d_ik| (over all genes and classes).
order_shrink <- function(delta, threshold) {
  dif <- rank(abs(delta)) - (length(delta) - threshold)
  dif <- matrix(dif, ncol = ncol(delta))
  delta <- delta * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}

shrink_fun <- function(method) {
  switch(method,
         soft  = soft_shrink,
         hard  = hard_shrink,
         order = order_shrink)
}

# ---- classification helpers -----------------------------------------------

# Row-wise arg-max of a discriminant-score matrix; ties go to the first class.
# (Named "softmax" in pamr.)
argmax_class <- function(x) {
  d <- dim(x)
  maxdist <- x[, 1]
  pclass <- rep(1, d[1])
  for (i in seq(2, d[2])) {
    l <- x[, i] > maxdist
    pclass[l] <- i
    maxdist[l] <- x[l, i]
  }
  dd <- dimnames(x)[[2]]
  if (is.null(dd) || !length(dd)) pclass
  else factor(pclass, levels = seq(d[2]), labels = dd)
}

safe_exp <- function(x) {
  xx <- sign(x) * pmin(abs(x), 500)
  exp(xx)
}

# Index of the minimum; among ties, the LAST one (i.e. the largest threshold
# for soft/hard thresholding, which gives the sparsest model).
last_which_min <- function(x) length(x) + 1 - which.min(rev(x))

# ---- cross-validation folds (from pamr) -----------------------------------

permute_rows <- function(x) {
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- stats::runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

#' Class-balanced cross-validation folds
#'
#' Splits sample indices into folds so that every class is spread as evenly as
#' possible across folds. This is the fold generator of the \pkg{pamr}
#' package, which is also used inside [hnsc_cv()].
#'
#' @param y Vector of class labels.
#' @param nfolds Number of folds. Defaults to the size of the smallest class,
#'   capped at 10 (and never less than 2).
#' @return A list of integer vectors of sample indices, one per fold.
#' @examples
#' set.seed(1)
#' balanced_folds(rep(c("A", "B"), c(6, 9)))
#' @export
balanced_folds <- function(y, nfolds = min(min(table(y)), 10)) {
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  nfolds <- max(nfolds, 2)
  yids <- split(seq(y), y)
  bigmat <- matrix(NA, ceiling(fmax / nfolds) * nfolds, length(totals))
  for (i in seq(totals)) {
    if (length(yids[[i]]) > 1) bigmat[seq(totals[i]), i] <- sample(yids[[i]])
    if (length(yids[[i]]) == 1) bigmat[seq(totals[i]), i] <- yids[[i]]
  }
  smallmat <- matrix(bigmat, nrow = nfolds)
  smallmat <- permute_rows(t(smallmat))
  res <- vector("list", nfolds)
  for (j in 1:nfolds) {
    jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j]
  }
  res
}

# ---- argument handling ------------------------------------------------------

method_table <- data.frame(
  label    = c("STh", "HTh", "OTh", "STh3", "HTh3", "OTh3"),
  method   = c("soft", "hard", "order", "soft", "hard", "order"),
  variance = c("homo", "homo", "homo", "hetero", "hetero", "hetero"),
  stringsAsFactors = FALSE
)

#' The six nearest shrunken centroid classifiers
#'
#' Lists the six classifiers implemented in the package, with the labels used
#' in the paper and the corresponding `method` / `variance` arguments.
#'
#' Only STh (soft thresholding with pooled variance) is the original PAM
#' classifier; HTh and OTh are the homoscedastic classifiers of Sahtout, Wang
#' and Ghimire (2024).
#'
#' @return A data frame with columns `label`, `method`, `variance`,
#'   `description` and `reference`.
#' @examples
#' hnsc_methods()
#' @export
hnsc_methods <- function() {
  out <- method_table
  out$description <- paste0(
    c("Soft", "Hard", "Order", "Soft", "Hard", "Order"), " thresholding, ",
    ifelse(out$variance == "homo", "pooled variance",
           "class-specific variances"))
  out$reference <- c("Tibshirani et al. (2002), original PAM",
                     "Sahtout, Wang and Ghimire (2024)",
                     "Sahtout, Wang and Ghimire (2024)",
                     "Sahtout and Wang", "Sahtout and Wang",
                     "Sahtout and Wang")
  out
}

method_label <- function(method, variance) {
  method_table$label[method_table$method == method &
                       method_table$variance == variance]
}

# Accepts either method = "soft"/"hard"/"order" together with variance, or a
# paper label such as "HTh3" in `method` (which then overrides `variance`).
resolve_method <- function(method, variance) {
  if (length(method) == 1 && method %in% method_table$label) {
    row <- method_table[method_table$label == method, ]
    return(list(method = row$method, variance = row$variance))
  }
  method <- match.arg(method[1], c("hard", "soft", "order"))
  variance <- match.arg(variance[1], c("hetero", "homo"))
  list(method = method, variance = variance)
}

check_xy <- function(x, y) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x))
    stop("`x` must be a numeric matrix with genes in rows and samples in columns.",
         call. = FALSE)
  if (length(y) != ncol(x)) {
    if (length(y) == nrow(x))
      stop("`x` must have genes in rows and samples in columns; ",
           "length(y) matches nrow(x), so try t(x).", call. = FALSE)
    stop("length(y) must equal ncol(x) (the number of samples).", call. = FALSE)
  }
  if (anyNA(x)) stop("`x` contains missing values.", call. = FALSE)
  if (anyNA(y)) stop("`y` contains missing values.", call. = FALSE)
  if (length(unique(y)) < 2) stop("`y` must contain at least two classes.",
                                  call. = FALSE)
  x
}
