make_data <- function(seed = 1) {
  set.seed(seed)
  simulate_hetero_data(ntr1 = 12, ntr2 = 15, nte1 = 20, nte2 = 20,
                       p = 200, a = 15)
}

test_that("simulation returns genes x samples matrices", {
  d <- make_data()
  expect_equal(dim(d$xtr), c(200, 27))
  expect_equal(dim(d$xte), c(200, 40))
  expect_equal(as.vector(table(d$ytr)), c(12, 15))
})

test_that("informative genes follow the paper's covariance", {
  set.seed(2)
  for (cc in c("B", "both")) {
    d <- simulate_hetero_data(400, 400, 1, 1, ngene = 5, p = 10,
                              s = rep(2, 10), a = 0, rho = 0.5,
                              correlated = cc)
    xb <- d$xtr[, d$ytr == "B"]
    xa <- d$xtr[, d$ytr == "A"]
    expect_gt(min(cor(t(xb[1:5, ]))), 0.999)           # rank-one Sigma
    expect_equal(unname(apply(xb[1:5, ], 1, var)), 4 * 0.5^(0:4),
                 tolerance = 0.2)                      # rho^(i-1) s_i^2
    expect_lt(abs(cor(xb[6, ], xb[7, ])), 0.2)          # other genes independent
    if (cc == "both") expect_gt(min(cor(t(xa[1:5, ]))), 0.999)
    else expect_lt(abs(cor(xa[1, ], xa[2, ])), 0.2)
  }
})

test_that("all six classifiers train, cross-validate and predict", {
  d <- make_data()
  for (m in hnsc_methods()$label) {
    fit <- hnsc_train(d$xtr, d$ytr, method = m)
    expect_s3_class(fit, "hnsc")
    expect_equal(fit$label, m)
    cv <- hnsc_cv(fit, d$xtr, d$ytr)
    expect_length(cv$error, length(fit$threshold))
    expect_true(all(cv$error >= 0 & cv$error <= 1))
    pr <- predict(fit, d$xte, threshold = fit$threshold[2])
    expect_length(pr, ncol(d$xte))
    expect_setequal(levels(pr), c("A", "B"))
    post <- predict(fit, d$xte, threshold = fit$threshold[2], type = "posterior")
    expect_equal(unname(rowSums(post)), rep(1, ncol(d$xte)))
  }
})

test_that("method labels and method/variance pairs are equivalent", {
  d <- make_data()
  a <- hnsc_train(d$xtr, d$ytr, method = "OTh3")
  b <- hnsc_train(d$xtr, d$ytr, method = "order", variance = "hetero")
  expect_equal(a$nonzero, b$nonzero)
  expect_equal(a$errors, b$errors)
})

test_that("hnsc selects a threshold and is reproducible with a seed", {
  d <- make_data()
  set.seed(5); f1 <- hnsc(d$xtr, d$ytr, method = "HTh3")
  set.seed(5); f2 <- hnsc(d$xtr, d$ytr, method = "HTh3")
  expect_equal(f1$selected.threshold, f2$selected.threshold)
  expect_equal(length(f1$selected.genes), f1$nonzero)
  expect_equal(predict(f1, d$xte), predict(f2, d$xte))
  cen <- predict(f1, type = "centroid")
  expect_equal(dim(cen), c(200, 2))
  f3 <- hnsc(d$xtr, d$ytr, method = "HTh3", threshold = 1)
  expect_null(f3$selection)
  expect_equal(f3$selected.threshold, 1)
})

test_that("input checks catch transposed data", {
  d <- make_data()
  expect_error(hnsc_train(t(d$xtr), d$ytr), "try t\\(x\\)")
  expect_error(hnsc_train(d$xtr, d$ytr[-1]), "length\\(y\\)")
})

test_that("heteroscedastic fit refuses singleton classes", {
  d <- make_data()
  y <- d$ytr; y[1] <- "C"
  expect_error(suppressWarnings(hnsc_train(d$xtr, y, variance = "hetero")),
               "at least 2")
})

test_that("metrics", {
  cm <- make_cm(c("A", "A", "B", "B"), c("A", "B", "B", "B"), c("A", "B"))
  expect_equal(unname(per_class_recall(cm)), c(0.5, 1))
  expect_equal(balanced_acc(cm), 0.75)
  d <- make_data()
  r <- compute_rho2(d$xtr, d$ytr)
  expect_true(is.logical(r$crossover_met))
  expect_equal(nrow(r), 2)
  expect_equal(r$threshold, round(2 / (c(12, 15) - 1), 4))
})

test_that("compare_nsc runs and summarises", {
  d <- make_data()
  res <- compare_nsc(d$xtr, d$ytr, d$xte, d$yte, n.runs = 1,
                     thresholding = "hard", verbose = FALSE)
  s <- summary(res)
  expect_equal(s$method, c("HTh", "HTh3"))
  expect_true(all(s$balanced.accuracy >= 0 & s$balanced.accuracy <= 1))
  rt <- recall_table(res)
  expect_equal(rt$class, c("A", "B"))
  expect_equal(mean(rt$recall_HTh), s$balanced.accuracy[1])
})
