data {
  int<lower = 0> N_obs;
  // location[N_obs] takes on value i for location number i etc.
  int location[N_obs];
  // days is how many days since prevalence exceeded a limit we have chosen
  vector[N_obs] days;
  // We model the daily number of newly diagnosed cases instead of cumulative numbers
  int new_cases[N_obs];
  
  int<lower = 0> N_locations;
  vector[N_locations] population;
}

transformed data {
  vector[N_obs] log_population = log(population[location]);
  vector[4] mu_priors;
  mu_priors[1] = 3.0;
  mu_priors[2] = -2.0;
  mu_priors[3] = 0.0;
  mu_priors[4] = -4.0;
}

parameters {
  
  cholesky_factor_corr[4] L_Omega;
  matrix[4, N_locations] pre_pars;
  vector<lower = 0>[4] sigma_pars;
  vector[4] mu_pars;
  
  real mu_phi_inv;
  real<lower = 0> sigma_phi_inv;
  vector[N_locations] log_phi_inv;
}

transformed parameters {
  vector[N_locations] log_alpha = to_vector(pre_pars[1]);
  vector<lower = 0>[N_locations] alpha = exp(log_alpha);
  
  
  vector[N_locations] log_beta = to_vector(pre_pars[2]);
  vector<lower = 0>[N_locations] beta = exp(log_beta);
  
  
  vector[N_locations] log_nu = to_vector(pre_pars[3]);
  vector<lower = 0>[N_locations] nu = exp(log_nu);
  
  
  vector[N_locations] logit_S = to_vector(pre_pars[4]);
  vector<lower = 0, upper = 1>[N_locations] S  = inv_logit(logit_S);
  vector<upper = 0>[N_locations] log_S = log(S);
  
  vector<lower = 0>[N_locations] phi_inv = exp(log_phi_inv);
  vector<lower = 0>[N_locations] phi = inv(phi_inv);
}

model {
  // Logistic equation calculations
  vector[N_obs] linear = beta[location] .* (days - alpha[location]);
  vector[N_obs] extra = (inv(nu[location]) + 1) .* log(1 + exp(log_nu[location] + linear));
  vector[N_obs] log_dfdt = log_beta[location] + log_S[location] + linear - extra;
  
  // Calculate covariance matrix for GLGM parameters
  matrix[4, 4] cov_mat = diag_pre_multiply(sigma_pars, L_Omega);
  
  // Multivariate Normal prior for the GLGM parameters
  mu_pars ~ normal(mu_priors, 1);
  sigma_pars ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2);
  
  for (i in 1:N_locations) pre_pars[, i] ~ multi_normal_cholesky(mu_pars, cov_mat);
  
  // Centered log-normal prior for the overdispersion
  mu_phi_inv ~ normal(-1, 1);
  sigma_phi_inv ~ exponential(1);
  log_phi_inv ~ normal(mu_phi_inv, sigma_phi_inv);
  
  //  Likelihood
  new_cases ~ neg_binomial_2_log(log_dfdt + log_population, phi[location]);
}

generated quantities {
  matrix[4, 4] corr_mat = multiply_lower_tri_self_transpose(L_Omega);
  real<lower = 0> sigma_alpha = sigma_pars[1];
  real<lower = 0> sigma_beta = sigma_pars[2];
  real<lower = 0> sigma_nu = sigma_pars[3];
  real<lower = 0> sigma_S = sigma_pars[4];
  real mu_alpha = mu_pars[1];
  real mu_beta = mu_pars[2];
  real mu_nu = mu_pars[3];
  real mu_S = mu_pars[4];
}
