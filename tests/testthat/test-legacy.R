# Checks that the package reproduces the original research code exactly.
# The original scripts are shipped in inst/legacy; the homoscedastic soft
# thresholding script needs pamr::pamr.predict.
test_that("results match the original research code", {
  skip_on_cran()
  skip_if_not_installed("pamr")
  leg <- system.file("legacy", package = "hetNSC")
  skip_if(leg == "")
  L <- new.env()
  for (f in c("pamr.hard.R", "pamr.soft.R", "pamr.order.R",
              "pamr.soft.hetero.R", "pamr.hard.hetero.R",
              "pamr.order.hetero.R"))
    sys.source(file.path(leg, f), envir = L)
  assign("pamr.predict", pamr::pamr.predict, envir = L)
  q <- function(expr) { utils::capture.output(r <- expr); r }

  set.seed(3)
  y <- rep(c("a", "b", "c"), c(8, 12, 10))
  gen <- function(lab) vapply(lab, function(l) {
    k <- match(l, c("a", "b", "c"))
    rnorm(150, mean = c(rep(k * 0.5, 10), rep(0, 140)), sd = 2 * k - 1)
  }, numeric(150))
  x <- gen(y)
  xt <- gen(rep(c("a", "b", "c"), 4))
  dat <- list(x = x, y = y)
  lg <- list(
    STh  = c("my.pamr.Scv.old", "pamr.train.old", "pamr.predict"),
    HTh  = c("my.pamr.Scv.hard", "pamr.train.hard", "pamr.predict.hard"),
    OTh  = c("my.pamr.Scv.order", "pamr.train.order", "pamr.predict.order"),
    STh3 = c("my.pamr.Scv.heter", "pamr.train.heter", "pamr.predict.heter"),
    HTh3 = c("my.pamr.Scv.H.heter", "pamr.train.H.heter", "pamr.predict.H.heter"),
    OTh3 = c("my.pamr.Scv.O.heter", "pamr.train.O.heter", "pamr.predict.O.heter"))
  for (m in names(lg)) {
    set.seed(77)
    pk <- q(L[[lg[[m]][1]]](dat))
    TR <- q(L[[lg[[m]][2]]](dat, threshold = pk[1]))
    post0 <- L[[lg[[m]][3]]](TR, xt, threshold = pk[1], type = "posterior")
    set.seed(77)
    fit <- hnsc(x, y, method = m)
    expect_equal(fit$selection$threshold, unname(pk[1]), info = m)
    expect_equal(fit$selection$error, unname(pk[3]), info = m)
    expect_equal(fit$nonzero, TR$nonzero, info = m)
    expect_equal(unname(predict(fit, xt, type = "posterior")), unname(post0),
                 info = m)
  }
})

test_that("simulated data match the original simulation scripts", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  leg <- system.file("legacy", package = "hetNSC")
  skip_if(leg == "")
  for (cfg in list(c("simu_heter_correlated.R", "B"),
                   c("simu_homo_correlated.R", "both"))) {
    L <- new.env()
    for (e in parse(file.path(leg, cfg[1])))
      if (is.call(e) && identical(e[[2]], as.name("generate_sample")))
        eval(e, L)
    environment(L$generate_sample) <- list2env(list(mvrnorm = MASS::mvrnorm),
                                               parent = globalenv())
    p <- 300
    s <- rep_len(1:20, p)
    set.seed(8900)
    d0 <- suppressMessages(L$generate_sample(10, 12, 5, 6, ngene = 20, p = p,
                                             s = s, a = 13))
    set.seed(8900)
    d1 <- simulate_hetero_data(10, 12, 5, 6, ngene = 20, p = p, s = s, a = 13,
                               correlated = cfg[2])
    expect_identical(unname(t(d0$Xtr)), unname(d1$xtr), info = cfg[1])
    expect_identical(unname(t(d0$Xte)), unname(d1$xte), info = cfg[1])
  }
})
