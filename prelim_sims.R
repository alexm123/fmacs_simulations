library(SimDesign)
library(lavaan)
library(pinsearch)
library(ggplot2)

# TODO:
#   - Increase the number of replications to 500
#   - Summarize the pattern of bias

# Define conditions: Testing different sample sizes
design <- createDesign(
  n = c(30, 100, 250)
)

# Fixed objects
set.seed(1855)

# Helper Function
# Generates a random covariance matrix
get_ucov <- function(p, scale = sqrt(0.1), n = 5) {
  W <- matrix(rnorm(p * n), nrow = n)
  WtW <- crossprod(W)
  D <- diag(1 / sqrt(diag(WtW))) * scale
  D %*% WtW %*% D
}
# Defines key parameters
fixed <- list(
  p = 6,
  lambda = c(.3, .7, .4, .5, .6, .4),
  dlambda = list(
    c(0, 0, 0, 0, 0, 0),
    c(.1, 0, 0, 0, 0, 0),
    c(.2, -.3, 0, 0, 0, 0),
    c(.3, -.3, 0, 0, 0, 0)
  ),
  nu = c(2, 3, 1.5, 3.5, 2, 3),
  alpha = c(0, -0.25, 0.25, 0.5),
  psi = c(1, 0.85, 1.15, 0.7),
  theta = c(1, 1.2, .8, .9, 1, 1) - .1,
  dtheta = matrix(
    runif(24, min = -0.2, max = 0.2),
    nrow = 4
  ),
  # ucov = replicate(4, get_ucov(6), simplify = FALSE)
  ucov = replicate(4, diag(.1, 6), simplify = FALSE),
  ninv_ind = c(1, 2)
)
# lavaan syntax for the measurement model
fixed$mod <- paste(
  "f =~",
  paste0("y", seq_len(fixed$p), collapse = " + ")
)
# Compute implied means and covariances
fixed <- within(fixed, {
  lambdag <- lapply(dlambda, FUN = \(x) x + lambda)
  Thetag <- lapply(seq_along(ucov),
                   FUN = function(g) {
                     diag(theta + dtheta[g, ]) + ucov[[g]]
                   })
  covy <- mapply(\(lam, psi, th) tcrossprod(lam) * psi + th,
                 lam = lambdag, psi = psi, th = Thetag,
                 SIMPLIFY = FALSE)
  meany <- mapply(\(lam, al, nu) nu + lam * al,
                  lam = lambdag, al = alpha, nu = list(nu),
                  SIMPLIFY = FALSE)
})

# Population effect size
fixed$fmacs_pop <- local({
  pooled_sd <- lapply(fixed$covy, FUN = \(x) diag(x)) |>
    do.call(what = rbind) |>
    colMeans() |>
    sqrt()

  fmacs(
    intercepts = matrix(rep(fixed$nu, 4),
                        nrow = 4,
                        byrow = TRUE
    ),
    loadings = sweep(
      do.call(rbind, fixed$dlambda),
      MARGIN = 2,
      STATS = fixed$lambda,
      FUN = "+"
    ),
    latent_mean = 0,
    latent_sd = 1,
    pooled_item_sd = pooled_sd
  )[1:2]
})

# Function for data generation
# sim_y <- function(n, lambda, nu, alpha, psi, Theta) {
#     covy <- tcrossprod(lambda) * psi + Theta
#     meany <- nu + lambda * alpha
#     MASS::mvrnorm(n, mu = meany, Sigma = covy)
# }
generate <- function(condition, fixed_objects) {
  ylist <- lapply(seq_along(fixed_objects$covy),
                  FUN = function(g) {
                    yg <- MASS::mvrnorm(
                      condition$n,
                      mu = fixed_objects$meany[[g]],
                      Sigma = fixed_objects$covy[[g]]
                    )
                    colnames(yg) <- paste0("y", seq_len(fixed_objects$p))
                    cbind(yg, group  = g)
                  })
  do.call(rbind, ylist)
}
sim1 <- generate(design[3, ], fixed_objects = fixed)


# Analysis
analyze <- function(condition, dat, fixed_objects) {
  # Define lavaan syntax
  pinv_fit <- cfa(
    fixed_objects$mod,
    data = dat,
    group = "group", std.lv = TRUE,
    group.equal = c("loadings", "intercepts"),
    group.partial = c(
      paste0("f=~y", fixed_objects$ninv_ind),
      paste0("y", fixed_objects$ninv_ind, "~1")
    )
  )
  as.vector(pinsearch::pin_effsize(pinv_fit))
}

# Evaluate/Summarize
evaluate <- function(condition, results, fixed_objects) {
  c(
    bias = colMeans(results) - fixed_objects$fmacs_pop,
    robust_bias = apply(results, 2, mean, trim = .1) -
      fixed_objects$fmacs_pop,
    emp_sd = apply(results, 2, sd),
    emp_mad = apply(results, 2, mad)
  )
}

out <- runSimulation(design,
                     replications = 500,
                     parallel = TRUE,
                     generate = generate,
                     analyse = analyze,
                     summarise = evaluate,
                     filename = "results-trial1",
                     packages = c("MASS", "lavaan", "pinsearch"),
                     fixed_objects = fixed
)


data <- readRDS("results-trial1.rds")
#write.csv(data, file = "out.csv")

# Bias for V1 and V2
ggplot(data, aes(x = n)) +
  geom_line(aes(y = bias.V1, color = "Bias V1")) +
  geom_point(aes(y = bias.V1, color = "Bias V1")) +
  geom_line(aes(y = bias.V2, color = "Bias V2")) +
  geom_point(aes(y = bias.V2, color = "Bias V2")) +
  labs(title = "Bias for V1 and V2 across Sample Sizes",
       x = "Sample Size", y = "Bias") +
  scale_color_manual(values = c("Bias V1" = "blue", "Bias V2" = "red")) +
  theme_minimal()

# Robust Bias for V1 and V2
ggplot(data, aes(x = n)) +
  geom_line(aes(y = robust_bias.V1, color = "Robust Bias V1")) +
  geom_point(aes(y = robust_bias.V1, color = "Robust Bias V1")) +
  geom_line(aes(y = robust_bias.V2, color = "Robust Bias V2")) +
  geom_point(aes(y = robust_bias.V2, color = "Robust Bias V2")) +
  labs(title = "Robust Bias for V1 and V2 across Sample Sizes",
       x = "Sample Size", y = "Robust Bias") +
  scale_color_manual(values = c("Robust Bias V1" = "blue", "Robust Bias V2" = "red")) +
  theme_minimal()

# Empirical Standard Deviation for V1 and V2
ggplot(data, aes(x = n)) +
  geom_line(aes(y = emp_sd.V1, color = "Empirical SD V1")) +
  geom_point(aes(y = emp_sd.V1, color = "Empirical SD V1")) +
  geom_line(aes(y = emp_sd.V2, color = "Empirical SD V2")) +
  geom_point(aes(y = emp_sd.V2, color = "Empirical SD V2")) +
  labs(title = "Empirical Standard Deviation for V1 and V2 across Sample Sizes",
       x = "Sample Size", y = "Empirical SD") +
  scale_color_manual(values = c("Empirical SD V1" = "blue", "Empirical SD V2" = "red")) +
  theme_minimal()

# Empirical MAD for V1 and V2
ggplot(data, aes(x = n)) +
  geom_line(aes(y = emp_mad.V1, color = "Empirical MAD V1")) +
  geom_point(aes(y = emp_mad.V1, color = "Empirical MAD V1")) +
  geom_line(aes(y = emp_mad.V2, color = "Empirical MAD V2")) +
  geom_point(aes(y = emp_mad.V2, color = "Empirical MAD V2")) +
  labs(title = "Empirical MAD for V1 and V2 across Sample Sizes",
       x = "Sample Size", y = "Empirical MAD") +
  scale_color_manual(values = c("Empirical MAD V1" = "blue", "Empirical MAD V2" = "red")) +
  theme_minimal()
