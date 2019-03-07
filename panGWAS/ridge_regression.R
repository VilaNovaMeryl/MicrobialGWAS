#' Helper function to compute the PC scores, estimates of SNP random effects, estimates of PC random effects and all
#' quantities required to perform a Wald test on PC random effects.
#'
#' @title pc_wald
#'
#' @param y Required. A named vector of phenotype. Two-level factors are automatically converted to binary responses.
#' @param x Required. A centered binary feature matrix. Typically produced by compressing a binary SNP genotype matrix
#' @param w Optional. Covariate matrix. Defaults to the intercept.
#' @param K Optional. Kinship matrix. Defaults to tcrossprod(x).
#' @param lambda_init Optional. Starting value for lambda (ratio between random effect variance and noise variance)
#' @param prefer.REML Optional. Logical. Should estimates of lambda and other quantities be computed using REML (TRUE, default) or MLE (FALSE) method.
#' @param skip.var Optional. Defaults to FALSE. Skip computation of posterior variance of SNP random effects?
#' @param maximize Optional. Maximize (TRUE, default) or evaluate (FALSE) the likelihood under lambda.
#'
#' @return A list with the following components \begin{itemize}
#' \item lambda_MLE Numeric. Estimate of lambda using the MLE method.
#' \item lambda_REML Numeric. Estimate of lambda using the REML method
#' \item a_MLE Numeric. Estimate of a (intercept) using the MLE method
#' \item a_REML Numeric. Estimate of a (intercept) using the REML method
#' \item tau_MLE Numeric. Estimate of tau (inverse variance of the noise) using the MLE method
#' \item tau_REML Numeric. Estimate of tau (inverse variance of the noise) using the REML method
#' \item Ebeta Numeric vector. Posterior mean of the SNP random effect
#' \item Vbeta Numeric vector. Posterior variance of the SNP random effect (only one value per SNP, diagonal of the full posterior variance covariance matrix)
#' \item Eg Numeric vector. Posterior mean of the PC random effect
#' \item Vg Numeric vector. Posterior variance of the PC random effect (stored as one value per PC, as the posterior variance covariance matrix is diagonal)
#' \item PC.scores Numeric matrix. PC scores of the strains, computed using a SVD of x
#' \item y.pred Numeric vector. Predicted phenotypes values (using posterior means of PC random effects)
#' @export
#'
#' @examples
pc_wald = function(y,                         ## Phenotype
                   x,                         ## Centered binary alleles matrix,
                   ##     can be reduced to speed up computations.
                   w = matrix(1,length(y),1), ## Default, intercept only
                   K = NULL,                  ## Kinship matrix
                   lambda_init = 200,
                   prefer.REML = TRUE,
                   skip.var = TRUE,           ## To speed up computations
                   maximize = TRUE            ## Maximize likelihood or simply evaluate it
) {
  if(lambda_init<=0) stop("ridge_regression: lambda_init must be positive")
  n = length(y)
  ## Sanity check on X
  if (any(abs(colMeans(x)) > sqrt(.Machine$double.eps))) stop("Matrix X should be centered.")
  X = x ## rename to fit with matrix algebra

  ## Convert factor/character discrete phenotype to integer
  if (is.character(y)) { y <- as.factor(y) }
  if (is.factor(y)) { y <- as.numeric(y) - 1 }

  # Precompute other quantities
  Wx = as.matrix(w)
  c = ncol(Wx)

  # Get Kinship matrix
  if (is.null(K)) {
    pattern.counts <- apply(X, 2, function(x) { diff(range(x))^2 })
    K <- tcrossprod(X) / sum(pattern.counts)
  }

  # Estimate lambda using maximum likelihood
  lambda_loglik = function(log_lambda) {
    lambda = exp(log_lambda)
    H = lambda * K + diag(n)
    Hinv = solve(H)
    Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
    ret = 0.5*n*log(0.5*n/pi) - 0.5*n - 0.5*as.numeric(determinant(H,log=TRUE)$modulus) - 0.5*n*log(t(y) %*% Px %*% y)
    if(!is.finite(ret)) return(-.Machine$double.xmax)
    return(ret)
  }
  if(maximize) {
    nlm_MLE = suppressWarnings(nlm(function(x) tryCatch(-lambda_loglik(x),error=function(e) -Inf),log(lambda_init)))
  } else {
    nlm_MLE = list("estimate" = log(lambda_init),
                   "minimum" = -lambda_loglik(log(lambda_init)))
  }
  lambda_MLE = exp(nlm_MLE$estimate)
  lambda_MLE <- min(lambda_MLE, 1e5) ## Hack
  ML = -nlm_MLE$minimum

  # Estimate lambda using REML
  lambda_logpartiallik = function(log_lambda) {
    lambda = exp(log_lambda)
    H = lambda * K + diag(n)
    Hinv = solve(H)
    Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
    # Note that whereas Zhou and Stephens use (n-c-1) we use (n-c) because we do not include an extra fixed effect for an individual genotype
    ret = 0.5*(n-c)*log(0.5*(n-c)/pi) - 0.5*(n-c) + 0.5*as.numeric(determinant(t(Wx)%*%Wx,log=TRUE)$modulus) - 0.5*as.numeric(determinant(H,log=TRUE)$modulus) - 0.5*as.numeric(determinant(t(Wx)%*%Hinv%*%Wx,log=TRUE)$modulus) - 0.5*(n-c)*log(t(y) %*% Px %*% y)
    if(!is.finite(ret)) return(-.Machine$double.xmax)
    return(ret)
  }
  if(maximize) {
    nlm_REML = suppressWarnings(nlm(function(x) tryCatch(-lambda_logpartiallik(x),error=function(e) -Inf),log(lambda_init)))
  } else {
    nlm_REML = list("estimate"=log(lambda_init),"minimum"=-lambda_logpartiallik(log(lambda_init)))
  }
  lambda_REML = exp(nlm_REML$estimate)
  lambda_REML <- min(lambda_REML, 1e5) ## Hack
  REML = -nlm_REML$minimum

  # Estimate a and tau given lambda using maximum likelihood
  lambda = lambda_MLE
  H = lambda * K + diag(n)
  Hinv = solve(H)
  a_MLE = c(solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv %*% y)
  Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
  tau_MLE = c(n / (t(y) %*% Px %*% y))

  # Estimate tau given lambda using REML
  lambda = lambda_REML
  H = lambda * K + diag(n)
  Hinv = solve(H)
  a_REML = c(solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv %*% y)
  Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
  # Note that whereas Zhou and Stephens use (n-c-1) we use (n-c) because we do not include an extra fixed effect for an individual genotype
  tau_REML = c((n - c) / (t(y) %*% Px %*% y))

  # Calculate the posterior mean and variance in effect size for each SNP
  lambda = lambda_REML
  tau = tau_REML
  a = a_REML
  if(!prefer.REML) {
    lambda = lambda_MLE
    tau = tau_MLE
    a = a_MLE
  }

  # Calculate posterior mean and variance of each snp mixed effect (gamma) in the ridge regression.

  # Note that with prior mean m=0, the posterior mean is the ridge regression estimate (Eq 11.57 of O`Hagan & Forster with tau=sigma^(-2))
  #	mstar	= solve(t(X) X + 1/lambda I_L) t(X) y
  # Using the svd decomposition of X = U D t(V) and the Woodburry formula
  # solve(lambda t(X) X + I_L) = (I_L - V [ lambda * D^2 (I_n + lambda * D^2)^{-1} ] t(V))
  # this simplifies to
  # mstar = V solve(D^2 + 1/lambda I_n) D t(U) y
  # TODO write corresponding matrix algebra
  svdX <- svd(X); U <- svdX$u; d <- svdX$d; V <- svdX$v
  ## Compute posterior mean of gamma
  post.mean <- c(V %*% (diag( lambda * d / (lambda * d^2 + 1) ) %*% t(U) %*% y))

  # Similarly, assuming a prior variance of W = lambda/tau I_L, the posterior variance is
  # Wstar = lambda/tau * (I_L + lambda t(X) X) = I_L - V [ lambda * D^2 (I_n + lambda * D^2)^{-1} ] t(V)
  # or Wstar = lambda/tau * trace(I_L - V %*% Cstar %*% t(V)) where
  # Cstar = (lambda * D^2 / (I_n + lambda * D^2))
  Cstar = diag(lambda * d^2 / (lambda * d^2 + 1))
  if(skip.var) {
    post.var = NA
  } else {
    post.var = lambda/tau * (1 - colSums(t(V^2)*diag(Cstar))) ## same as lambda/tau * diag(I_L - V %*% Cstar %*% t(V)), we don't compute off-diagonal elements of the posterior matrix.
  }

  # Calculate mean and posterior variance of the PC mixed effect (g) using the MLE estimates for lambda and tau
  # After PCA, T = X %*% rotation with rotation = V and T = UD from svd decomposition X = UD t(V)
  # Note that the PC mixed effect is simply g = t(rotation) %*% gamma = t(V) %*% gamma
  # gstar = t(rotation) %*% mstar = solve(D^2 + 1/lambda I_n) D t(U) y
  # Sstar = t(rotation) %*% Wstar %*% rotation = lambda/tau * (I_n - Cstar) = lambda/tau * (I_n + lambda * D^2)^{-1}
  pca.post.mean <- c(diag( lambda * d / (lambda * d^2 + 1)) %*% t(U) %*% y) ## n x 1
  pca.post.var  <- (lambda / tau) * 1 / (lambda * d^2 + 1) ## take only diagonal as Sstar is a diagonal matrix
  ## Wald test: compare pca.post.mean^2 / pca.post.var to chisq with 1 df

  ## Predicted phenotype values using effect sizes
  # \hat{y} = a + X \hat{u} = a + T \hat{g} = U D gstar with
  # - \hat{g} = solve(D^2 + 1/lambda I_n) D t(U) y
  # - T = U D
  y.pred <- a + U %*% diag( lambda * d^2 / (lambda * d^2 + 1)) %*% t(U) %*% y # n x 1

  ## PC.scores
  PC.scores <- U %*% diag(d)
  rownames(PC.scores) <- rownames(X)
  colnames(PC.scores) <- paste0("PC", 1:ncol(PC.scores))

  ret = list("lambda_MLE"  = lambda_MLE,
             "lambda_REML" = lambda_REML,
             "a_MLE"       = a_MLE,
             "a_REML"      = a_REML,
             "tau_MLE"     = tau_MLE,
             "tau_REML"    = tau_REML,
             "Ebeta"       = post.mean,
             "Vbeta"       = post.var,
             "Eg"          = pca.post.mean,
             "Vg"          = pca.post.var,
             "PC.scores"   = PC.scores,
             "y.pred"      = y.pred
             )
  return(ret)
}
