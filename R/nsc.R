# Core computational engine shared by all six classifiers.
#
# Homoscedastic branch  = pamr's nsc(): soft shrinkage is the original PAM (STh);
#                         hard / order shrinkage are HTh, OTh (Sahtout, Wang and Ghimire 2024).
# Heteroscedastic branch = nsc.heter(), nsc.H.heter(), nsc.O.heter() of
#                         Sahtout and Wang (STh3, HTh3, OTh3).
#
# The arithmetic is kept in the same order as the original code so that the
# results agree with it to machine precision.

# Discriminant scores for the pooled-variance classifier (pamr's diag.disc).
# x and centroids are already standardised by (s_i + s_0).
# Returns -delta_k(x*)/2 up to a class-independent constant (larger = better).
diag_disc <- function(x, centroids, prior, weight) {
  if (!missing(weight)) {
    posid <- (weight > 0)
    if (any(posid)) {
      weight <- sqrt(weight[posid])
      centroids <- centroids[posid, , drop = FALSE] * weight
      x <- x[posid, , drop = FALSE] * weight
    } else {
      mat <- outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  dd <- t(x) %*% centroids
  dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2)) / 2 - log(prior)
  names(dd0) <- NULL
  scale(dd, dd0, FALSE)
}

# Discriminant scores for the heteroscedastic classifier, eq. (2.9) of the
# paper. x = x* - xbar (unscaled); centroids = shrunken centroid offsets
# xbar'_k - xbar; sik2 = s_ik^2 + s_0k^2 (p x K).
# Returns -delta_k(x*)/2 (larger = better).
diag_disc_heter <- function(x, centroids, prior, weight, sik2) {
  if (!missing(weight)) {
    posid <- (weight > 0)
    if (any(posid)) {
      weight <- sqrt(weight[posid])
      centroids <- centroids[posid, , drop = FALSE] * weight
      x <- x[posid, , drop = FALSE] * weight
      sik2 <- sik2[posid, , drop = FALSE] * (weight^2)
    } else {
      mat <- outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  centroids2.scaled <- centroids^2 / sik2
  dd0 <- drop(rep(1, nrow(centroids)) %*% centroids2.scaled) / 2 - log(prior) +
    drop(rep(1, nrow(sik2)) %*% log(sqrt(2 * pi * sik2)))
  names(dd0) <- NULL
  dd <- t(x) %*% (centroids / sik2)
  K <- ncol(centroids)
  disc <- NULL
  for (k in 1:K) {
    dd3 <- t((x^2) / sik2[, k]) %*% rep(1, nrow(x)) / 2
    disc <- cbind(disc, dd3)
  }
  scale(dd - disc, dd0, FALSE)
}

# Fit an NSC classifier over a grid of thresholds.
#
# x, xtest: p x n matrices.  y: factor.  If xtest is NULL the training data
# are also used as test data (resubstitution error).  se.scale may be passed
# in (cross-validation re-uses the full-data value, as pamr does).
nsc_fit <- function(x, y, xtest = NULL, ytest = NULL,
                    method = "hard", variance = "hetero",
                    threshold = NULL, n.threshold = 30,
                    offset.percent = 50, prior = NULL, se.scale = NULL,
                    remove.zeros = TRUE, sign.contrast = "both",
                    compute.prob = TRUE) {
  if (!is.factor(y)) y <- factor(y)
  argy <- if (is.null(ytest)) y else ytest
  n.class <- table(y)
  if (is.null(prior)) {
    prior <- table(y) / length(y)
    prior <- prior[prior != 0]
  }
  if (min(n.class) == 1)
    warning("a class contains only 1 sample", call. = FALSE)
  if (variance == "hetero" && min(n.class) < 2)
    stop("The heteroscedastic classifiers need at least 2 training samples ",
         "per class.", call. = FALSE)
  if (is.null(xtest)) {
    xtest <- x
    ytest <- y
  }
  n <- sum(n.class)
  ntest <- ncol(xtest)
  K <- length(prior)
  p <- nrow(x)
  Y <- stats::model.matrix(~ factor(y) - 1, data = list(y = y))
  dimnames(Y) <- list(NULL, names(n.class))

  # class centroids  xbar_ik  (p x K)
  centroids <- scale(x %*% Y, FALSE, n.class)
  xdif <- x - centroids %*% t(Y)
  threshold.scale <- rep(1, K)
  names(threshold.scale) <- names(n.class)
  sik2 <- NULL

  if (variance == "homo") {
    # pooled within-class sd  s_i + s_0   (eq. 2.3)
    sd <- (xdif^2) %*% rep(1 / (n - K), n)
    sd <- drop(sqrt(sd))
    offset <- stats::quantile(sd, offset.percent / 100)
    sd <- sd + offset
    if (is.null(se.scale)) se.scale <- sqrt(1 / n.class - 1 / n)
  } else {
    # class-specific variances s_ik^2 + s_0k^2   (eq. 2.6-2.7)
    sik2 <- scale((xdif^2) %*% Y, FALSE, n.class - 1)
    offset <- apply(sik2, 2, function(v) stats::quantile(v, offset.percent / 100))
    sik2 <- scale(sik2, -offset, FALSE)
    # denominator of d_ik, eq. (2.7)
    sd <- scale(sik2, FALSE, 1 / (1 / n.class - 2 / n)) +
      drop(scale(sik2, FALSE, 1 / n.class) %*% rep(1 / n^2, K))
    sd <- sqrt(sd)
    if (is.null(se.scale)) se.scale <- 1
  }

  centroid.overall <- drop(x %*% rep(1 / n, n))
  delta <- (centroids - centroid.overall) / sd
  if (variance == "homo")
    delta <- scale(delta, FALSE, threshold.scale * se.scale)
  if (sign.contrast == "positive") delta <- delta * (delta > 0)
  if (sign.contrast == "negative") delta <- delta * (delta < 0)

  if (!is.null(threshold)) {
    n.threshold <- length(threshold)
  } else if (method == "order") {
    mx <- round(length(delta) / 4)
    threshold <- round(seq(mx, 1, length = n.threshold))
  } else {
    threshold <- seq(0, max(abs(delta)), length = n.threshold)
  }

  shrink <- shrink_fun(method)
  nonzero <- seq(n.threshold)
  errors <- threshold
  yhat <- matrix(NA_character_, ntest, n.threshold)
  prob <- if (compute.prob) array(0, c(ntest, K, n.threshold)) else NULL
  xs <- if (variance == "homo") (xtest - centroid.overall) / sd
        else xtest - centroid.overall

  for (ii in 1:n.threshold) {
    delta.shrunk <- shrink(delta, threshold[ii])
    nonzero[ii] <- attr(delta.shrunk, "nonzero")
    if (variance == "homo") {
      delta.shrunk <- scale(delta.shrunk, FALSE, 1 / (threshold.scale * se.scale))
      posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
      dd <- diag_disc(xs, delta.shrunk, prior, weight = posid)
    } else {
      posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
      dd <- diag_disc_heter(xs, delta.shrunk * sd, prior, weight = posid, sik2)
    }
    cls <- argmax_class(dd)
    yhat[, ii] <- as.character(cls)
    if (compute.prob) {
      dd <- safe_exp(dd)
      prob[, , ii] <- dd / drop(dd %*% rep(1, K))
    }
    if (!is.null(ytest)) errors[ii] <- sum(cls != ytest)
  }

  if (remove.zeros) n.threshold <- match(0, nonzero, n.threshold)
  keep <- seq(n.threshold)
  if (compute.prob) {
    dimnames(prob) <- list(colnames(xtest), names(n.class),
                           format(round(threshold, 3)))
    prob <- prob[, , keep, drop = FALSE]
  }
  object <- list(
    y = argy, yhat = yhat[, keep, drop = FALSE], prob = prob,
    centroids = centroids, centroid.overall = centroid.overall,
    sd = sd, sik2 = sik2, offset = offset, delta = delta,
    threshold = threshold[keep], nonzero = nonzero[keep],
    threshold.scale = threshold.scale, se.scale = se.scale,
    prior = prior, sign.contrast = sign.contrast,
    method = method, variance = variance,
    label = method_label(method, variance),
    offset.percent = offset.percent
  )
  if (!is.null(ytest)) object$errors <- errors[keep]
  object
}

# Shrunken quantities of a fitted object at one threshold.
shrunken_delta <- function(fit, threshold) {
  K <- length(fit$prior)
  if (fit$variance == "homo") {
    delta <- scale((fit$centroids - fit$centroid.overall) / fit$sd, FALSE,
                   fit$threshold.scale * fit$se.scale)
  } else {
    delta <- (fit$centroids - fit$centroid.overall) / fit$sd
  }
  if (fit$sign.contrast == "positive") delta <- delta * (delta > 0)
  if (fit$sign.contrast == "negative") delta <- delta * (delta < 0)
  ds <- shrink_fun(fit$method)(delta, threshold)
  if (fit$variance == "homo")
    ds <- scale(ds, FALSE, 1 / (fit$threshold.scale * fit$se.scale))
  ds
}

# Discriminant scores (log-posterior up to a constant), n_new x K.
nsc_scores <- function(fit, newx, threshold) {
  ds <- shrunken_delta(fit, threshold)
  K <- length(fit$prior)
  posid <- drop(abs(ds) %*% rep(1, K)) > 0
  if (fit$variance == "homo") {
    diag_disc((newx - fit$centroid.overall) / fit$sd, ds, fit$prior, posid)
  } else {
    diag_disc_heter(newx - fit$centroid.overall, ds * fit$sd, fit$prior,
                    weight = posid, fit$sik2)
  }
}
