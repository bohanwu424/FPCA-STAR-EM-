STAR_FPCA_EM = function (y, estimator_em = james_estimator,estimator_initialize = estimator.fpca, transformation = "np", y_max = Inf, 
          sd_init = 10, tol = 10^-10, max_iters = 1000,show_iter = FALSE,smooth_Lam = 5000) 
{
  if (any(y < 0) || any(y != floor(y))) 
    stop("y must be nonnegative counts")
  if (any(y > y_max)) 
    stop("y must not exceed y_max")
  transformation = tolower(transformation)
  if (!is.element(transformation, c("identity", "log", "sqrt", 
                                    "np", "pois", "neg-bin", "box-cox"))) 
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', 'neg-bin', or 'box-cox'")
  transform_family = ifelse(test = is.element(transformation, 
                                              c("identity", "log", "sqrt", "box-cox")), yes = "bc", 
                            no = "cdf")
  n = length(y)
  if (transform_family == "bc") {
    if (transformation == "identity") 
      lambda = 1
    if (transformation == "log") 
      lambda = 0
    if (transformation == "sqrt") 
      lambda = 1/2
    if (transformation == "box-cox") 
      lambda = runif(n = 1)
    g = function(t) g_bc(t, lambda = lambda)
    g_inv = function(s) g_inv_bc(s, lambda = lambda)
    sum_log_deriv = (lambda - 1) * sum(log(y + 1))
  }
  if (transform_family == "cdf") {
    g = g_cdf(y = y, distribution = transformation)
    t_grid = sort(unique(round(c(seq(0, min(2 * max(y), 
                                            y_max), length.out = 250), quantile(unique(y[y < 
                                                                                           y_max] + 1), seq(0, 1, length.out = 250))), 8)))
    g_inv = g_inv_approx(g = g, t_grid = t_grid)
    sum_log_deriv = sum(log(pmax(g(y + 1, deriv = 1), 0.01)))
    lambda = NULL
  }
  z_hat = g(y + 1)
  fit = estimator_initialize(z_hat)#intialize
  if (is.null(fit$fitted.values) || is.null(fit$coefficients)) 
    stop("The estimator() function must return 'fitted.values' and 'coefficients'")
  mu_hat = fit$fitted.values
  theta_hat = fit$coefficients
  sigma_hat = sd(z_hat - mu_hat)
  logLik0 = logLik_em0 = sum_log_deriv + sum(dnorm(z_hat, 
                                                   mean = mu_hat, sd = sigma_hat, log = TRUE))
  if (sd_init > 0) {
    z_hat = g(y + 1) + sd_init * sigma_hat * rnorm(n = n)
    fit = estimator_initialize(z_hat)
    mu_hat = fit$fitted.values
    theta_hat = fit$coefficients
    sigma_hat = sd(z_hat - mu_hat)
  }
  p = length(theta_hat)
  a_y = a_j(y, y_max = y_max)
  a_yp1 = a_j(y + 1, y_max = y_max)
  z_lower = g(a_y)
  z_upper = g(a_yp1)
  mu_all = zhat_all = array(0, c(max_iters, n))
  theta_all = array(0, c(max_iters, p))
  sigma_all = numeric(max_iters)
  logLik_all = numeric(max_iters)
  iter_all = numeric(max_iters)
  for (s in 1:max_iters) {
    z_mom = truncnorm_mom(a = z_lower, b = z_upper, mu = mu_hat, 
                          sig = sigma_hat)
    z_hat = z_mom$m1
    z2_hat = z_mom$m2
    if (any(is.infinite(z_hat))) {
      warning("Infinite z_hat values: returning the problematic indices")
      return(list(error_inds = which(is.infinite(z_hat))))
    }
    fit = estimator_em(z_hat,theta_hat,smooth_Lam)
    mu_hat = fit$fitted.values
    theta_hat = fit$coefficients
    if(show_iter){
      print(fit$iter)
    }
    sigma_hat = sqrt((sum(z2_hat) + sum(mu_hat^2) - 2 * 
                        sum(z_hat * mu_hat))/n)
    if (transformation == "box-cox") {
      ff <- function(l_bc) {
        sapply(l_bc, function(l_bc) {
          -logLikeRcpp(g_a_j = g_bc(a_y, lambda = l_bc), 
                       g_a_jp1 = g_bc(a_yp1, lambda = l_bc), mu = mu_hat, 
                       sigma = rep(sigma_hat, n))
        })
      }
      a = 0
      b = 1
      while (ff(b) == Inf) {
        b = b * 0.8
      }
      lambda = BrentMethod(a, b, fcn = ff, tol = .Machine$double.eps^0.2)$x
      g = function(t) g_bc(t, lambda = lambda)
      g_inv = function(s) g_inv_bc(s, lambda = lambda)
      z_lower = g(a_y)
      z_upper = g(a_yp1)
    }
    logLik_em = logLikeRcpp(g_a_j = z_lower, g_a_jp1 = z_upper, 
                            mu = mu_hat, sigma = rep(sigma_hat, n))
    mu_all[s, ] = mu_hat
    theta_all[s, ] = theta_hat
    sigma_all[s] = sigma_hat
    logLik_all[s] = logLik_em
    zhat_all[s, ] = z_hat
    iter_all[s] = fit$iter
    if ((logLik_em - logLik_em0)^2 < tol) 
      break
    logLik_em0 = logLik_em
  }
  mu_all = mu_all[1:s, ]
  theta_all = theta_all[1:s, ]
  sigma_all = sigma_all[1:s]
  logLik_all = logLik_all[1:s]
  zhat_all = zhat_all[1:s, ]
  iter_all = iter_all[1:s]
  if (y_max < Inf) {
    Jmax = rep(y_max + 1, n)
  }
  else {
    Jmax = round_fun(g_inv(qnorm(0.9999, mean = mu_hat, 
                                 sd = sigma_hat)), y_max = y_max)
    Jmax[Jmax > 2 * max(y)] = 2 * max(y)
  }
  Jmaxmax = max(Jmax)
  y_hat = expectation_gRcpp(g_a_j = g(a_j(0:Jmaxmax, y_max = y_max)), 
                            g_a_jp1 = g(a_j(1:(Jmaxmax + 1), y_max = y_max)), mu = mu_hat, 
                            sigma = rep(sigma_hat, n), Jmax = Jmax)
  resids_ds = qnorm(runif(n) * (pnorm((z_upper - mu_hat)/sigma_hat) - 
                                  pnorm((z_lower - mu_hat)/sigma_hat)) + pnorm((z_lower - 
                                                                                  mu_hat)/sigma_hat))
  resids_ds_rep = sapply(1:10, function(...) qnorm(runif(n) * 
                                                     (pnorm((z_upper - mu_hat)/sigma_hat) - pnorm((z_lower - 
                                                                                                     mu_hat)/sigma_hat)) + pnorm((z_lower - mu_hat)/sigma_hat)))
  list(coefficients = theta_hat, fitted.values = y_hat, g.hat = g, 
       sigma.hat = sigma_hat, mu.hat = mu_hat, z.hat = z_hat, 
       residuals = resids_ds, residuals_rep = resids_ds_rep, 
       logLik = logLik_em, logLik0 = logLik0, lambda = lambda, 
       mu_all = mu_all, theta_all = theta_all, sigma_all = sigma_all, 
       logLik_all = logLik_all, zhat_all = zhat_all, transformation = transformation, 
       y_max = y_max, tol = tol, max_iters = max_iters, iter_all = iter_all)
}
