library(here)
library(tidyverse)
library(knitr)
library(kableExtra)
library(broom)
library(cowplot)
library(rstan)
library(tidybayes)
library(scales)
library(lubridate)
library(plotly)
library(posterior)

theme_set(theme_classic(base_size = 12) + 
              background_grid(minor = "none", major = "none") +
              theme(legend.position = "none"))

fitted_dates <- here("Past Models", "Data") %>% 
    list.files %>% 
    str_match("_(2020-[0-9]{2}-[0-9]{2})\\.csv") %>% 
    .[, 2] %>% 
    na.omit %>% 
    ymd

read_models <- function(read_date) {
    here("Past Models", "Models", str_c("Model_", read_date, ".rds")) %>% 
        read_rds()
}

read_data <- function(read_date) {
    here("Past Models", "Data", str_c("COVID_Data_", read_date, ".csv")) %>% 
        read_csv(col_types = cols())
}

convergence_info <- function(model) {
    model %>% 
        as_draws %>% 
        as_draws_df %>% 
        as_tibble %>% 
        pivot_longer(c(-.chain, -.iteration, -.draw)) %>% 
        select(-.iteration, -.draw) %>% 
        group_by(.chain, name) %>% 
        mutate(iter = row_number()) %>% 
        ungroup %>% 
        pivot_wider(names_from = .chain, values_from = value) %>% 
        rename(chain = 3:6) %>% 
        group_by(name) %>% 
        summarise(bulk = ess_bulk(cbind(chain1, chain2, chain3, chain4)),
                  tail = ess_tail(cbind(chain1, chain2, chain3, chain4)),
                  ess_mean = ess_mean(cbind(chain1, chain2, chain3, chain4)),
                  ess_sd = ess_sd(cbind(chain1, chain2, chain3, chain4)),
                  rhat = rhat(cbind(chain1, chain2, chain3, chain4))) %>% 
        summarise(bulk = min(bulk, na.rm = TRUE),
                  tail = min(tail, na.rm = TRUE),
                  ess_mean = min(ess_mean, na.rm = TRUE),
                  ess_sd = min(ess_sd, na.rm = TRUE),
                  rhat = max(rhat, na.rm = TRUE))
}

m <- map(fitted_dates, read_models)

d <- map(fitted_dates, read_data)

conv_info <- m %>% 
    map_dfr(convergence_info)

tab <- m %>% 
    map(get_elapsed_time) %>% 
    map(as_tibble, rownames = "chain") %>% 
    reduce(bind_rows) %>% 
    mutate(fit_date = fitted_dates[cumsum(chain == "chain:1")],
           total_time = warmup + sample) %>% 
    group_by(fit_date) %>% 
    summarise(max_time = max(total_time))

n_div <- m %>%  map_dbl(get_num_divergent)

n_leapfrog <- m %>% map_dbl(~ get_num_leapfrog_per_iteration(.) %>% mean %>% round)

n_info <- d %>% 
    map_dfr(~ summarise(., n_obs = n(), n_loc = length(unique(location))))


tab <- tab %>% 
    mutate(fit_date = format(fit_date, "%B %d")) %>% 
    bind_cols(n_info) %>% 
    bind_cols(conv_info) %>% 
    mutate(n_div = n_div,
           n_leapfrog = n_leapfrog) %>% 
    mutate(max_time = round(max_time))

tab
