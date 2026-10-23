#' Train a nearest shrunken centroid classifier over a grid of thresholds
#'
#' Fits one of the six nearest shrunken centroid (NSC) classifiers for a
#' sequence of threshold values. `variance = "homo"` uses the PAM test
#' statistic and discriminant function with a pooled within-class standard
#' deviation; combined with `method = "soft"` this is the original PAM
#' classifier (STh, Tibshirani et al. 2002), while `method = "hard"` or
#' `"order"` gives the homoscedastic classifiers of Sahtout, Wang and Ghimire
#' (2024) (HTh, OTh). `variance = "hetero"` gives the heteroscedastic statistic
#' of Sahtout and Wang with class-specific variances in both the test
#' statistic (eq. 2.7) and the discriminant function (eq. 2.9) (STh3, HTh3,
#' OTh3).
#'
#' Most users will want [hnsc()], which also chooses the threshold by
#' cross-validation.
#'
#' @param x Numeric matrix of expression values with **genes in rows and
#'   samples in columns** (the same layout as \pkg{pamr}).
#' @param y Vector or factor of class labels, one per column of `x`.
#' @param method Thresholding rule: `"hard"`, `"soft"` or `"order"`. A paper
#'   label (`"STh"`, `"HTh"`, `"OTh"`, `"STh3"`, `"HTh3"`, `"OTh3"`) is also
#'   accepted, in which case `variance` is ignored. See [hnsc_methods()].
#' @param variance `"hetero"` (class-specific variances, the proposed method)
#'   or `"homo"` (pooled variance). With `"homo"`, `method = "soft"` is the
#'   original PAM classifier; `"hard"` and `"order"` are the homoscedastic
#'   classifiers of Sahtout, Wang and Ghimire (2024).
#' @param threshold Optional vector of thresholds. For soft and hard
#'   thresholding these are cut-offs on |d_ik|; for order thresholding they are
#'   the number of largest |d_ik| values (over all genes and classes) that are
#'   kept. If `NULL`, a default grid of `n.threshold` values is used.
#' @param n.threshold Number of thresholds in the default grid.
#' @param offset.percent Percentile of the standard deviations (homoscedastic)
#'   or of the class-specific variances (heteroscedastic) added as the
#'   stabilising constant s_0 (s_0k). Default 50, the median.
#' @param prior Class prior probabilities. Defaults to the class proportions
#'   in the training data.
#' @param sign.contrast `"both"` (default), `"positive"` or `"negative"`:
#'   restrict to centroid differences of one sign.
#' @param remove.zeros Drop thresholds beyond the first one that keeps no genes.
#'
#' @return An object of class `"hnsc"`: a list containing, among others,
#'   `threshold`, `nonzero` (number of genes kept at each threshold),
#'   `errors` (training resubstitution errors), `centroids`,
#'   `centroid.overall`, `sd` (denominator of d_ik), `sik2`
#'   (s_ik^2 + s_0k^2, heteroscedastic only), `delta` (the unshrunken test
#'   statistics d_ik), `prior`, `method`, `variance` and `label`.
#' @seealso [hnsc()], [hnsc_cv()], [predict.hnsc()]
#' @examples
#' set.seed(1)
#' dat <- simulate_hetero_data(ntr1 = 10, ntr2 = 10, nte1 = 5, nte2 = 5,
#'                             p = 200, a = 15)
#' fit <- hnsc_train(dat$xtr, dat$ytr, method = "hard", variance = "hetero")
#' fit
#' @export
hnsc_train <- function(x, y, method = c("hard", "soft", "order"),
                       variance = c("hetero", "homo"),
                       threshold = NULL, n.threshold = 30,
                       offset.percent = 50, prior = NULL,
                       sign.contrast = c("both", "positive", "negative"),
                       remove.zeros = TRUE) {
  x <- check_xy(x, y)
  m <- resolve_method(method, variance)
  sign.contrast <- match.arg(sign.contrast)
  y <- factor(y)
  if (is.null(prior)) {
    prior <- table(y) / length(y)
    prior <- prior[prior != 0]
  }
  fit <- nsc_fit(x, y, method = m$method, variance = m$variance,
                 threshold = threshold, n.threshold = n.threshold,
                 offset.percent = offset.percent, prior = prior,
                 remove.zeros = remove.zeros, sign.contrast = sign.contrast)
  fit$genenames <- rownames(x)
  fit$call <- match.call()
  class(fit) <- "hnsc"
  fit
}

#' Cross-validation for a fitted NSC classifier
#'
#' Estimates the misclassification error at each threshold of a fitted
#' [hnsc_train()] object by class-balanced K-fold cross-validation (the
#' procedure of \pkg{pamr}'s `pamr.cv`). The number of folds is the size of
#' the smallest class, capped at 10.
#'
#' @param fit An object returned by [hnsc_train()].
#' @param x,y The training data used to create `fit`.
#' @param folds Optional list of sample-index vectors defining the folds
#'   (for example from [balanced_folds()]).
#' @return A list of class `"hnsc_cv"` with components `threshold`, `size`
#'   (number of genes kept by the full-data fit), `error` (cross-validated
#'   error rate), `yhat` (matrix of cross-validated predictions) and `folds`.
#' @examples
#' set.seed(1)
#' dat <- simulate_hetero_data(ntr1 = 10, ntr2 = 10, nte1 = 5, nte2 = 5,
#'                             p = 200, a = 15)
#' fit <- hnsc_train(dat$xtr, dat$ytr, method = "HTh3")
#' cv <- hnsc_cv(fit, dat$xtr, dat$ytr)
#' cv
#' @export
hnsc_cv <- function(fit, x, y, folds = NULL) {
  x <- check_xy(x, y)
  cv_core(fit, x, factor(y), folds = folds)
}

cv_core <- function(fit, x, y, folds = NULL) {
  n <- length(y)
  if (is.null(folds)) folds <- balanced_folds(y)
  threshold <- fit$threshold
  yhat <- matrix(as.character(y), n, length(threshold))
  for (ii in seq_along(folds)) {
    f <- folds[[ii]]
    a <- nsc_fit(x[, -f, drop = FALSE], y[-f], x[, f, drop = FALSE],
                 method = fit$method, variance = fit$variance,
                 threshold = threshold, offset.percent = fit$offset.percent,
                 prior = fit$prior, se.scale = fit$se.scale,
                 remove.zeros = FALSE, sign.contrast = fit$sign.contrast,
                 compute.prob = FALSE)
    yhat[f, ] <- a$yhat
  }
  error <- colSums(yhat != as.character(y)) / n
  out <- list(threshold = threshold, size = fit$nonzero, error = error,
              yhat = yhat, y = y, folds = folds, label = fit$label)
  class(out) <- "hnsc_cv"
  out
}

#' Choose the threshold by refined ("deep") cross-validation
#'
#' Implements the threshold search of Sahtout, Wang and Ghimire (2024) used in
#' the paper: a cross-validation over the default threshold grid, followed by
#' repeated refinement of the grid between the neighbours of the current best
#' threshold until the number of selected genes can no longer change. Among
#' thresholds with equal cross-validated error, the one examined last (for
#' soft and hard thresholding typically the larger, sparser one) is chosen.
#'
#' @inheritParams hnsc_train
#' @param max.iter Safety cap on the number of refinement rounds.
#' @return A list of class `"hnsc_select"` with `threshold` (selected value),
#'   `size` (genes kept), `error` (its CV error), `path` (matrix of every
#'   threshold / size / CV error evaluated), `method`, `variance` and `label`.
#' @examples
#' set.seed(1)
#' dat <- simulate_hetero_data(ntr1 = 10, ntr2 = 10, nte1 = 5, nte2 = 5,
#'                             p = 200, a = 15)
#' sel <- hnsc_select(dat$xtr, dat$ytr, method = "OTh3")
#' sel$threshold
#' @export
hnsc_select <- function(x, y, method = c("hard", "soft", "order"),
                        variance = c("hetero", "homo"),
                        offset.percent = 50, max.iter = 100) {
  x <- check_xy(x, y)
  m <- resolve_method(method, variance)
  deep_cv(x, factor(y), m$method, m$variance, offset.percent, max.iter)
}

deep_cv <- function(x, y, method, variance, offset.percent = 50,
                    max.iter = 100) {
  prior <- table(y) / length(y)
  prior <- prior[prior != 0]
  train <- function(threshold = NULL)
    nsc_fit(x, y, method = method, variance = variance, threshold = threshold,
            offset.percent = offset.percent, prior = prior)

  fit0 <- train()
  cv0 <- cv_core(fit0, x, y)
  path <- cbind(threshold = cv0$threshold, size = cv0$size, error = cv0$error)
  nr <- nrow(path)
  loca <- last_which_min(cv0$error)
  lo <- max(1, loca - 1)
  hi <- min(loca + 1, nr)
  nonGenDif <- cv0$size[lo] - cv0$size[hi]
  if (nonGenDif > 1) {
    q1 <- seq(cv0$threshold[lo], cv0$threshold[hi],
              length.out = min(30, nonGenDif))
    iter <- 0
    repeat {
      iter <- iter + 1
      fit <- train(q1)
      cv <- cv_core(fit, x, y)
      sm <- cbind(threshold = cv$threshold, size = cv$size, error = cv$error)
      path <- rbind(path, sm)
      where <- last_which_min(sm[, 3])
      lo <- max(1, where - 1)
      hi <- min(where + 1, nrow(sm))
      cond1 <- sm[lo, 2] - sm[hi, 2]
      if (!(cond1 > 1 && nrow(sm) > 2) || iter >= max.iter) break
      q1 <- seq(sm[lo, 1], sm[hi, 1], length.out = min(30, cond1))
    }
  }
  best <- path[last_which_min(path[, 3]), ]
  out <- list(threshold = unname(best[1]), size = unname(best[2]),
              error = unname(best[3]), path = path,
              method = method, variance = variance,
              label = method_label(method, variance))
  class(out) <- "hnsc_select"
  out
}

#' Fit a nearest shrunken centroid classifier with the threshold chosen by
#' cross-validation
#'
#' The main user-facing function. It selects the threshold with the refined
#' cross-validation search of [hnsc_select()] (unless `threshold` is given)
#' and then fits the classifier on all training data at that threshold.
#' Use [predict.hnsc()] to classify new samples.
#'
#' @inheritParams hnsc_train
#' @param threshold Optional single threshold. If `NULL` (default) it is
#'   chosen by [hnsc_select()].
#' @param max.iter Passed to [hnsc_select()].
#' @return An object of class `"hnsc"` (see [hnsc_train()]) fitted at the
#'   selected threshold, with an extra component `selection` holding the
#'   [hnsc_select()] result (or `NULL` if `threshold` was supplied) and
#'   `selected.genes`, the row indices of the genes used by the classifier.
#' @references
#' Sahtout, M. O. and Wang, H. Heteroscedastic nearest shrunken centroid
#' classification: methods and applications to high-dimensional data.
#' *Statistical Applications in Genetics and Molecular Biology* (accepted).
#'
#' Sahtout, M. O., Wang, H. and Ghimire, S. (2024). Different thresholding
#' methods on nearest shrunken centroid algorithm. *Communications in
#' Statistics - Simulation and Computation*, 53(3), 1444-1460.
#' \doi{10.1080/03610918.2022.2047201}
#'
#' Tibshirani, R., Hastie, T., Narasimhan, B. and Chu, G. (2002). Diagnosis of
#' multiple cancer types by shrunken centroids of gene expression. *PNAS*,
#' 99(10), 6567-6572.
#' @examples
#' set.seed(1)
#' dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 50, nte2 = 50,
#'                             p = 500, a = 15)
#' fit <- hnsc(dat$xtr, dat$ytr, method = "HTh3")
#' fit
#' pred <- predict(fit, dat$xte)
#' mean(pred != dat$yte)            # test error
#' head(predict(fit, dat$xte, type = "posterior"))
#' @export
hnsc <- function(x, y, method = c("hard", "soft", "order"),
                 variance = c("hetero", "homo"), threshold = NULL,
                 offset.percent = 50, prior = NULL, max.iter = 100) {
  x <- check_xy(x, y)
  m <- resolve_method(method, variance)
  sel <- NULL
  if (is.null(threshold)) {
    sel <- deep_cv(x, factor(y), m$method, m$variance, offset.percent, max.iter)
    threshold <- sel$threshold
  }
  if (length(threshold) != 1) stop("`threshold` must be a single value.",
                                   call. = FALSE)
  fit <- hnsc_train(x, y, method = m$method, variance = m$variance,
                    threshold = threshold, offset.percent = offset.percent,
                    prior = prior, remove.zeros = FALSE)
  fit$selection <- sel
  fit$selected.threshold <- threshold
  fit$selected.genes <- predict.hnsc(fit, type = "nonzero")
  fit$call <- match.call()
  fit
}

#' Predict with a nearest shrunken centroid classifier
#'
#' @param object An object returned by [hnsc()] or [hnsc_train()].
#' @param newx Matrix of new samples, genes in rows (same genes and order as
#'   the training data). Not needed for `type = "centroid"` or `"nonzero"`.
#' @param threshold Threshold to use. Defaults to the threshold selected by
#'   [hnsc()]; must be supplied for objects from [hnsc_train()].
#' @param type `"class"` (predicted labels), `"posterior"` (class posterior
#'   probabilities), `"centroid"` (shrunken class centroids, genes x classes)
#'   or `"nonzero"` (row indices of the genes used).
#' @param ... Not used.
#' @return A factor, matrix or integer vector depending on `type`.
#' @examples
#' set.seed(1)
#' dat <- simulate_hetero_data(ntr1 = 10, ntr2 = 10, nte1 = 5, nte2 = 5,
#'                             p = 200, a = 15)
#' fit <- hnsc_train(dat$xtr, dat$ytr, method = "STh3")
#' predict(fit, dat$xte, threshold = fit$threshold[5])
#' @export
predict.hnsc <- function(object, newx = NULL, threshold = NULL,
                         type = c("class", "posterior", "centroid", "nonzero"),
                         ...) {
  type <- match.arg(type)
  if (is.null(threshold)) threshold <- object$selected.threshold
  if (is.null(threshold))
    stop("Supply `threshold` (this object was not fitted by hnsc()).",
         call. = FALSE)
  if (length(threshold) != 1) stop("`threshold` must be a single value.",
                                   call. = FALSE)
  if (type %in% c("class", "posterior")) {
    if (is.null(newx)) stop("`newx` is required for type = '", type, "'.",
                            call. = FALSE)
    if (is.data.frame(newx)) newx <- as.matrix(newx)
    if (is.vector(newx)) newx <- matrix(newx, ncol = 1)
    if (nrow(newx) != nrow(object$centroids))
      stop("`newx` must have the same genes (rows) as the training data.",
           call. = FALSE)
    dd <- nsc_scores(object, newx, threshold)
    if (type == "class") {
      cls <- argmax_class(dd)
      names(cls) <- colnames(newx)
      return(cls)
    }
    dd <- safe_exp(dd)
    post <- dd / drop(dd %*% rep(1, length(object$prior)))
    rownames(post) <- colnames(newx)
    return(post)
  }
  ds <- shrunken_delta(object, threshold)
  if (type == "centroid") {
    cen <- object$centroid.overall + ds * object$sd
    attributes(cen) <- list(dim = dim(cen),
                            dimnames = list(object$genenames,
                                            colnames(object$centroids)))
    return(cen)
  }
  nz <- drop(abs(ds) %*% rep(1, ncol(object$centroids))) > 0
  seq(nz)[nz]
}

#' @export
print.hnsc <- function(x, ...) {
  cat("Nearest shrunken centroid classifier:", x$label,
      sprintf("(%s thresholding, %s variance)\n", x$method,
              if (x$variance == "hetero") "class-specific" else "pooled"))
  cat("Classes:", paste(colnames(x$centroids), collapse = ", "),
      " | genes:", nrow(x$centroids), "\n")
  if (!is.null(x$selected.threshold)) {
    cat(sprintf("Selected threshold: %s  (genes used: %d",
                format(signif(x$selected.threshold, 4)),
                length(x$selected.genes)))
    if (!is.null(x$selection))
      cat(sprintf(", CV error: %.4f", x$selection$error))
    cat(")\n")
  } else {
    tab <- data.frame(threshold = signif(x$threshold, 4),
                      nonzero = x$nonzero,
                      train.errors = x$errors)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.hnsc_cv <- function(x, ...) {
  cat("Cross-validation for", x$label, "\n")
  print(data.frame(threshold = signif(x$threshold, 4), nonzero = x$size,
                   cv.error = round(x$error, 4)), row.names = FALSE)
  invisible(x)
}

#' @export
print.hnsc_select <- function(x, ...) {
  cat("Deep CV threshold search for", x$label, "\n")
  cat(sprintf("Selected threshold: %s | genes: %d | CV error: %.4f",
              format(signif(x$threshold, 4)), as.integer(x$size), x$error),
      sprintf("(%d thresholds evaluated)\n", nrow(x$path)))
  invisible(x)
}
