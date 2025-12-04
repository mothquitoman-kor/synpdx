# R/models.R
# Internal model utilities

#' @keywords internal
.safe_term <- function(k, t) ifelse(abs(k) < 1e-8, t, (1 - exp(-k * t)) / k)

#' @keywords internal
.gomp_y <- function(alpha1, alpha2, t) {
  ifelse(abs(alpha2) < 1e-8, alpha1 * t, (alpha1 / alpha2) * (1 - exp(-alpha2 * t)))
}

#' @keywords internal
.logistic_y <- function(K, r, t) {
  log(K) - log(1 + (K - 1) * exp(-r * t))
}

# -------------------------------------------------------------------
# Exponential growth + asymmetric drug effects
# -------------------------------------------------------------------

# exp + const A, const B
#' @keywords internal
exp_constA_constB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  a      <- exp(psi[id, 1])
  beta_a <- exp(psi[id, 2])
  beta_b <- exp(psi[id, 3])
  a * t - (beta_a * IA + beta_b * IB) * t
}

# exp + const A, decay B
#' @keywords internal
exp_constA_decayB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  a      <- exp(psi[id, 1])
  beta_a <- exp(psi[id, 2])
  beta_b <- exp(psi[id, 3])
  k_b    <- exp(psi[id, 4])
  a * t - beta_a * IA * t - beta_b * IB * .safe_term(k_b, t)
}

# exp + decay A, const B
#' @keywords internal
exp_decayA_constB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  a      <- exp(psi[id, 1])
  beta_a <- exp(psi[id, 2])
  beta_b <- exp(psi[id, 3])
  k_a    <- exp(psi[id, 4])
  a * t - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * t
}

# exp + decay A, decay B
#' @keywords internal
exp_decayA_decayB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  a      <- exp(psi[id, 1])
  beta_a <- exp(psi[id, 2])
  beta_b <- exp(psi[id, 3])
  k_a    <- exp(psi[id, 4])
  k_b    <- exp(psi[id, 5])
  a * t - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * .safe_term(k_b, t)
}

# -------------------------------------------------------------------
# Gompertz growth + asymmetric drug effects
# -------------------------------------------------------------------

# gomp + const A, const B
#' @keywords internal
gomp_constA_constB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  alpha1 <- exp(psi[id, 1])
  alpha2 <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  .gomp_y(alpha1, alpha2, t) - (beta_a * IA + beta_b * IB) * t
}

# gomp + const A, decay B
#' @keywords internal
gomp_constA_decayB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  alpha1 <- exp(psi[id, 1])
  alpha2 <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  k_b    <- exp(psi[id, 5])
  .gomp_y(alpha1, alpha2, t) - beta_a * IA * t - beta_b * IB * .safe_term(k_b, t)
}

# gomp + decay A, const B
#' @keywords internal
gomp_decayA_constB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  alpha1 <- exp(psi[id, 1])
  alpha2 <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  k_a    <- exp(psi[id, 5])
  .gomp_y(alpha1, alpha2, t) - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * t
}

# gomp + decay A, decay B
#' @keywords internal
gomp_decayA_decayB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  alpha1 <- exp(psi[id, 1])
  alpha2 <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  k_a    <- exp(psi[id, 5])
  k_b    <- exp(psi[id, 6])
  .gomp_y(alpha1, alpha2, t) - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * .safe_term(k_b, t)
}

# -------------------------------------------------------------------
# Logistic growth + asymmetric drug effects
# -------------------------------------------------------------------

# logistic + const A, const B
#' @keywords internal
logistic_constA_constB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  K      <- exp(psi[id, 1])
  r      <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  .logistic_y(K, r, t) - (beta_a * IA + beta_b * IB) * t
}

# logistic + const A, decay B
#' @keywords internal
logistic_constA_decayB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  K      <- exp(psi[id, 1])
  r      <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  k_b    <- exp(psi[id, 5])
  .logistic_y(K, r, t) - beta_a * IA * t - beta_b * IB * .safe_term(k_b, t)
}

# logistic + decay A, const B
#' @keywords internal
logistic_decayA_constB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  K      <- exp(psi[id, 1])
  r      <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  k_a    <- exp(psi[id, 5])
  .logistic_y(K, r, t) - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * t
}

# logistic + decay A, decay B
#' @keywords internal
logistic_decayA_decayB <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  K      <- exp(psi[id, 1])
  r      <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3])
  beta_b <- exp(psi[id, 4])
  k_a    <- exp(psi[id, 5])
  k_b    <- exp(psi[id, 6])
  .logistic_y(K, r, t) - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * .safe_term(k_b, t)
}

# -------------------------------------------------------------------
# Initial values & parameter names
# -------------------------------------------------------------------

attr(exp_constA_constB,      "theta0")   <- c(log(0.05), log(0.2), log(0.2))
attr(exp_constA_constB,      "parnames") <- c("a","beta_a","beta_b")

attr(exp_constA_decayB,      "theta0")   <- c(log(0.05), log(0.2), log(0.2), log(0.1))
attr(exp_constA_decayB,      "parnames") <- c("a","beta_a","beta_b","k_b")

attr(exp_decayA_constB,      "theta0")   <- c(log(0.05), log(0.2), log(0.2), log(0.1))
attr(exp_decayA_constB,      "parnames") <- c("a","beta_a","beta_b","k_a")

attr(exp_decayA_decayB,      "theta0")   <- c(log(0.05), log(0.2), log(0.2), log(0.1), log(0.1))
attr(exp_decayA_decayB,      "parnames") <- c("a","beta_a","beta_b","k_a","k_b")

attr(gomp_constA_constB,     "theta0")   <- c(log(0.1), log(0.02), log(0.2), log(0.2))
attr(gomp_constA_constB,     "parnames") <- c("alpha1","alpha2","beta_a","beta_b")

attr(gomp_constA_decayB,     "theta0")   <- c(log(0.1), log(0.02), log(0.2), log(0.2), log(0.1))
attr(gomp_constA_decayB,     "parnames") <- c("alpha1","alpha2","beta_a","beta_b","k_b")

attr(gomp_decayA_constB,     "theta0")   <- c(log(0.1), log(0.02), log(0.2), log(0.2), log(0.1))
attr(gomp_decayA_constB,     "parnames") <- c("alpha1","alpha2","beta_a","beta_b","k_a")

attr(gomp_decayA_decayB,     "theta0")   <- c(log(0.1), log(0.02), log(0.2), log(0.2), log(0.1), log(0.1))
attr(gomp_decayA_decayB,     "parnames") <- c("alpha1","alpha2","beta_a","beta_b","k_a","k_b")

attr(logistic_constA_constB, "theta0")   <- c(log(10), log(0.1), log(0.2), log(0.2))
attr(logistic_constA_constB, "parnames") <- c("K","r","beta_a","beta_b")

attr(logistic_constA_decayB, "theta0")   <- c(log(10), log(0.1), log(0.2), log(0.2), log(0.1))
attr(logistic_constA_decayB, "parnames") <- c("K","r","beta_a","beta_b","k_b")

attr(logistic_decayA_constB, "theta0")   <- c(log(10), log(0.1), log(0.2), log(0.2), log(0.1))
attr(logistic_decayA_constB, "parnames") <- c("K","r","beta_a","beta_b","k_a")

attr(logistic_decayA_decayB, "theta0")   <- c(log(10), log(0.1), log(0.2), log(0.2), log(0.1), log(0.1))
attr(logistic_decayA_decayB, "parnames") <- c("K","r","beta_a","beta_b","k_a","k_b")
