#' Simulate two-class data with correlated informative genes
#'
#' Generates the two-class designs with correlated informative genes used in
#' the paper's simulation study. Class A has mean zero for every gene; class B
#' has mean `a/10 * i` for the first `ngene` (informative) genes and zero
#' elsewhere. All non-informative genes are independent with standard
#' deviation `s[i]` in both classes.
#'
#' The informative genes are drawn jointly from a multivariate normal
#' distribution with covariance
#' \deqn{\Sigma = v v^\top, \qquad v_i = \sqrt{\rho^{\,i-1}}\, s_i,}{Sigma = v v', v_i = sqrt(rho^(i-1)) s_i,}
#' i.e. \eqn{\Sigma_{ij} = \rho^{(i+j-2)/2} s_i s_j}{Sigma_ij = rho^((i+j-2)/2) s_i s_j}.
#' This is the matrix used to produce the results in the paper. It has rank
#' one, so the informative genes are perfectly correlated with one another,
#' and gene `i` has variance \eqn{\rho^{i-1} s_i^2}{rho^(i-1) s_i^2}.
#'
#' @param ntr1,ntr2 Training sample sizes of classes A and B.
#' @param nte1,nte2 Test sample sizes of classes A and B.
#' @param ngene Number of informative genes.
#' @param p Total number of genes.
#' @param s Vector of length `p` with the gene standard deviations. The
#'   default cycles through 1, 2, ..., 20 as in the paper.
#' @param a Signal strength; class-B means are `a/10 * (1:ngene)`.
#' @param rho Correlation parameter (0.95 in the paper).
#' @param correlated Which classes have correlated informative genes:
#'   * `"B"` (default): only class B; class A is independent with
#'     standard deviation `s[i]` for every gene. This is the heteroscedastic
#'     correlated design (script `simu_heter_correlated.R`).
#'   * `"both"`: both classes share the covariance above, so the design is
#'     homoscedastic (script `simu_homo_correlated.R`).
#' @return A list with `xtr`, `xte` (genes x samples matrices), `ytr`, `yte`
#'   (character label vectors with values `"A"` and `"B"`).
#' @details Random numbers are drawn in the same order as in the original
#'   scripts, so with the same seed the data are identical (after
#'   transposition) to those used for the paper.
#' @examples
#' set.seed(8900)
#' dat <- simulate_hetero_data(ntr1 = 15, ntr2 = 20, nte1 = 50, nte2 = 50,
#'                             p = 1000, a = 12)
#' dim(dat$xtr)
#' table(dat$ytr)
#' @export
simulate_hetero_data <- function(ntr1, ntr2, nte1, nte2, ngene = 20,
                                 p = 10000,
                                 s = rep_len(1:20, p),
                                 a = 10, rho = 0.95,
                                 correlated = c("B", "both")) {
  correlated <- match.arg(correlated)
  stopifnot(length(s) == p, ngene >= 1, ngene < p)
  N1 <- ntr1 + nte1
  N2 <- ntr2 + nte2
  inf <- seq_len(ngene)
  other <- (ngene + 1):p

  mu2 <- a / 10 * inf
  cors <- rho^(0:(ngene - 1))
  Sigma <- matrix(sqrt(cors) * s[inf], ncol = 1) %*%
    matrix(sqrt(cors) * s[inf], nrow = 1)

  if (correlated == "B") {
    # order of random draws as in simu_heter_correlated.R
    x2_other <- scale(matrix(stats::rnorm(N2 * (p - ngene), 0, 1),
                             nrow = N2, ncol = p - ngene),
                      center = FALSE, scale = 1 / s[other])
    x2_inf <- MASS::mvrnorm(n = N2, mu = mu2, Sigma = Sigma)
    x1 <- scale(matrix(stats::rnorm(p * N1, mean = 0, sd = 1),
                       nrow = N1, ncol = p),
                center = FALSE, scale = 1 / s)
  } else {
    # order of random draws as in simu_homo_correlated.R
    x1_inf <- MASS::mvrnorm(n = N1, mu = rep(0, ngene), Sigma = Sigma)
    x2_inf <- MASS::mvrnorm(n = N2, mu = mu2, Sigma = Sigma)
    x1_other <- scale(matrix(stats::rnorm((p - ngene) * N1, mean = 0, sd = 1),
                             nrow = N1, ncol = p - ngene),
                      center = FALSE, scale = 1 / s[other])
    x2_other <- scale(matrix(stats::rnorm(N2 * (p - ngene), 0, 1),
                             nrow = N2, ncol = p - ngene),
                      center = FALSE, scale = 1 / s[other])
    x1 <- cbind(matrix(x1_inf, nrow = N1), x1_other)
  }
  x2 <- cbind(matrix(x2_inf, nrow = N2), x2_other)
  attributes(x1) <- list(dim = dim(x1))
  attributes(x2) <- list(dim = dim(x2))

  Xtr <- rbind(x1[1:ntr1, , drop = FALSE], x2[1:ntr2, , drop = FALSE])
  Xte <- rbind(x1[-(1:ntr1), , drop = FALSE], x2[-(1:ntr2), , drop = FALSE])
  genes <- paste0("g", seq_len(p))
  xtr <- t(Xtr)
  xte <- t(Xte)
  dimnames(xtr) <- list(genes, paste0("tr", seq_len(ncol(xtr))))
  dimnames(xte) <- list(genes, paste0("te", seq_len(ncol(xte))))
  list(xtr = xtr, ytr = c(rep("A", ntr1), rep("B", ntr2)),
       xte = xte, yte = c(rep("A", nte1), rep("B", nte2)))
}
