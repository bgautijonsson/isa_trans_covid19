data {
  int<lower = 0> N_obs;
  // location[N_obs] takes on value i for location number i etc.
  int location[N_obs];
  // days is how many days since case rate per 1000 exceeded a limit we have chosen
  vector[N_obs] days;
  // We model the daily number of newly diagnosed cases instead of cumulative numbers
  int new_cases[N_obs];
  
  int<lower = 0> N_locations;
  vector[N_locations] population;
}

transformed data {
  vector[N_obs] log_population = log(population[location]);
  vector[N_locations] log_population_cent = log(population) - mean(log(population));
  vector[4] mu_priors;
  mu_priors[1] = 3.0;
  mu_priors[2] = -2.0;
  mu_priors[3] = 0.0;
  mu_priors[4] = -4.0;
}

parameters {
  
  cholesky_factor_corr[4] L_Omega;
  matrix[4, N_locations] z_pars;
  vector<lower = 0>[4] sigma_pars;
  real mu_alpha;
  real mu_beta;
  real mu_nu;
  real mu_S;
  
  real mu_phi_inv;
  real<lower = 0> sigma_phi_inv;
  vector[N_locations] z_phi_inv;
}

transformed parameters {
  matrix[4, N_locations] cov_mat = diag_pre_multiply(sigma_pars, L_Omega) * z_pars;
  
  vector[N_locations] log_alpha = to_vector(cov_mat[1]) + mu_alpha;
  vector[N_locations] alpha = exp(log_alpha);
  
  
  vector[N_locations] log_beta = to_vector(cov_mat[2]) + mu_beta;
  vector[N_locations] beta = exp(log_beta);
  
  
  vector[N_locations] log_nu = to_vector(cov_mat[3]) + mu_nu;
  vector[N_locations] nu = exp(log_nu);
  

  vector[N_locations] logit_S = to_vector(cov_mat[4]) + mu_S;
  vector<lower = 0, upper = 1>[N_locations] S = inv_logit(logit_S);
  vector<upper = 0>[N_locations] log_S = log(S);
  
  
  vector[N_locations] log_phi_inv = mu_phi_inv + sigma_phi_inv * z_phi_inv;
  vector<lower = 0>[N_locations] phi_inv = exp(log_phi_inv);
  vector<lower = 0>[N_locations] phi = inv(phi_inv);
}

model {
  // Logistic equation calculations
  vector[N_obs] linear = beta[location] .* (days - alpha[location]);
  vector[N_obs] extra = (inv(nu[location]) + 1) .* log(1 + exp(log_nu[location] + linear));
  vector[N_obs] log_dfdt = log_beta[location] + log_S[location] + linear - extra;
  
  // Multivariate Normal prior for the GLGM parameters
  mu_alpha ~ normal(3, 1);
  mu_beta ~ normal(-2, 1);
  mu_nu ~ normal(0, 1);
  mu_S ~ normal(-4, 1);
  sigma_pars ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z_pars) ~ std_normal();
  
  // Non-Centered prior for the overdispersion
  mu_phi_inv ~ normal(-1, 1);
  sigma_phi_inv ~ exponential(1);
  z_phi_inv ~ std_normal();
  
  //  Likelihood
  new_cases ~ neg_binomial_2_log(log_dfdt + log_population, phi[location]);
}

generated quantities {
  matrix[4, 4] corr_mat = multiply_lower_tri_self_transpose(L_Omega);
  real<lower = 0> sigma_alpha = sigma_pars[1];
  real<lower = 0> sigma_beta = sigma_pars[2];
  real<lower = 0> sigma_nu = sigma_pars[3];
  real<lower = 0> sigma_S = sigma_pars[4];
  // Logistic equation calculations
  // vector[N_obs] linear = beta[location] .* (days - alpha[location]);
  // vector[N_obs] extra = (exp(-log_nu[location]) + 1) .* log(1 + exp(log_nu[location] + linear));
  // vector[N_obs] log_dfdt = log_beta[location] + log_S[location] + linear - extra;
  // vector[N_locations] log_beta_tilde = log_beta - log_nu;
  // vector<lower = 0>[N_locations] beta_tilde = exp(log_beta_tilde);
  // real log_lik[N_obs];
  // int sim_cases[N_obs]= neg_binomial_2_rng(dfdt .* population[location], phi[location]);
  // for (i in 1:N_obs) log_lik[i] = neg_binomial_2_log_lpmf(new_cases[i] | log_dfdt[i] + log_population[i], phi[location[i]]);
}
