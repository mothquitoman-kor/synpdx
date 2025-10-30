# models.R
#' Internal model utilities
#' @keywords internal
.safe_term <- function(k, t) ifelse(abs(k) < 1e-8, t, (1 - exp(-k * t)) / k)

#' @keywords internal
.gomp_y <- function(alpha1, alpha2, t) ifelse(abs(alpha2) < 1e-8, alpha1 * t, (alpha1/alpha2) * (1 - exp(-alpha2 * t)))

#' @keywords internal
.logistic_y <- function(K, r, t) log(K) - log(1 + (K - 1) * exp(-r * t))

# Exponential growth + constant effects
#' @keywords internal
exp_const <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  a      <- exp(psi[id, 1]); beta_a <- exp(psi[id, 2]); beta_b <- exp(psi[id, 3])
  a * t - (beta_a * IA + beta_b * IB) * t
}

# Exponential growth + decaying effects
#' @keywords internal
exp_decay <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  a      <- exp(psi[id, 1]); beta_a <- exp(psi[id, 2]); beta_b <- exp(psi[id, 3])
  k_a    <- exp(psi[id, 4]); k_b    <- exp(psi[id, 5])
  a * t - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * .safe_term(k_b, t)
}

# Gompertz + constant
#' @keywords internal
gomp_const <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  alpha1 <- exp(psi[id, 1]); alpha2 <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3]); beta_b <- exp(psi[id, 4])
  .gomp_y(alpha1, alpha2, t) - (beta_a * IA + beta_b * IB) * t
}

# Gompertz + decay
#' @keywords internal
gomp_decay <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  alpha1 <- exp(psi[id, 1]); alpha2 <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3]); beta_b <- exp(psi[id, 4])
  k_a    <- exp(psi[id, 5]); k_b    <- exp(psi[id, 6])
  .gomp_y(alpha1, alpha2, t) - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * .safe_term(k_b, t)
}

# Logistic + constant
#' @keywords internal
logistic_const <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  K <- exp(psi[id, 1]); r <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3]); beta_b <- exp(psi[id, 4])
  .logistic_y(K, r, t) - (beta_a * IA + beta_b * IB) * t
}

# Logistic + decay
#' @keywords internal
logistic_decay <- function(psi, id, xidep) {
  t  <- xidep$time; IA <- xidep$I_A; IB <- xidep$I_B
  K <- exp(psi[id, 1]); r <- exp(psi[id, 2])
  beta_a <- exp(psi[id, 3]); beta_b <- exp(psi[id, 4])
  k_a    <- exp(psi[id, 5]); k_b    <- exp(psi[id, 6])
  .logistic_y(K, r, t) - beta_a * IA * .safe_term(k_a, t) - beta_b * IB * .safe_term(k_b, t)
}


# --- initial values & names ---
attr(exp_const,      "theta0")   <- c(log(0.05), log(0.2),  log(0.2))
attr(exp_const,      "parnames") <- c("a","beta_a","beta_b")

attr(exp_decay,      "theta0")   <- c(log(0.05), log(0.2),  log(0.2),  log(0.1), log(0.1))
attr(exp_decay,      "parnames") <- c("a","beta_a","beta_b","k_a","k_b")

attr(gomp_const,     "theta0")   <- c(log(0.1),  log(0.02), log(0.2),  log(0.2))
attr(gomp_const,     "parnames") <- c("alpha1","alpha2","beta_a","beta_b")

attr(gomp_decay,     "theta0")   <- c(log(0.1),  log(0.02), log(0.2),  log(0.2),  log(0.1), log(0.1))
attr(gomp_decay,     "parnames") <- c("alpha1","alpha2","beta_a","beta_b","k_a","k_b")

attr(logistic_const, "theta0")   <- c(log(10),   log(0.1),  log(0.2),  log(0.2))
attr(logistic_const, "parnames") <- c("K","r","beta_a","beta_b")

attr(logistic_decay, "theta0")   <- c(log(10),   log(0.1),  log(0.2),  log(0.2),  log(0.1), log(0.1))
attr(logistic_decay, "parnames") <- c("K","r","beta_a","beta_b","k_a","k_b")
