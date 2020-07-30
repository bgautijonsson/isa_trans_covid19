library(here)
library(tidyverse)
library(lubridate)
library(tidybayes)

fitted_dates <- here("Past Models", "Models") %>% 
    list.files %>% 
    str_match("_(2020-[0-9]{2}-[0-9]{2})\\.rds") %>% 
    .[, 2] %>% 
    na.omit %>% 
    ymd


for (i in seq_along(fitted_dates)) {
    
    if (file.exists(here("Past Models", "Predictions", str_c("Predictions_", fitted_dates[i], ".csv")))) {
        next
    }
    
    writeLines(str_c("\n",
                     "Predicting for model from: ", 
                     fitted_dates[i],
                     "\n"))
    
    
    m <- here("Past Models", "Models", str_c("Model_", fitted_dates[i], ".rds")) %>% 
        read_rds
    
    d <- here("Past Models", "Data", str_c("COVID_Data_", fitted_dates[i], ".csv")) %>% 
        read_csv
    
    res <- spread_draws(m, 
                        alpha[location_id],
                        beta[location_id], 
                        S[location_id],
                        phi[location_id],
                        nu[location_id], 
                        n = 1000) %>% 
        ungroup %>% 
        select(-.chain, .iteration, .draw) %>% 
        group_by(location_id) %>% 
        mutate(iter = row_number()) %>% 
        ungroup %>% 
        inner_join(
            d %>% 
                group_by(location, location_id, population) %>% 
                summarise(start_date = min(date),
                          start_cases = min(total_cases)),
            by = "location_id"
        ) %>% 
        expand_grid(date = seq.Date(ymd("2020-02-01"), ymd("2020-09-01"), by = 1)) %>% 
        filter(date > start_date) %>% 
        mutate(days = as.numeric(date - start_date)) %>% 
        select(-start_date) %>% 
        mutate(z = beta * (days - alpha),
               f = S - S / (1 + nu * exp(z))^(1/nu),
               dfdt = beta * (S - f) * (1 - ((S - f) / S)^nu) / nu,
               new_cases = rnbinom(n(), mu = dfdt * population, size = phi)) %>% 
        group_by(location, iter) %>% 
        mutate(total_cases = as.numeric(cumsum(new_cases * (days > 0))) + start_cases) %>% 
        ungroup %>% 
        pivot_longer(c(new_cases, total_cases), names_to = "variable", values_to = "value") %>% 
        group_by(location, location_id, date, variable) %>% 
        summarise(median = median(value),
                  conf.low = quantile(value, .025),
                  conf.high = quantile(value, .975))
    
    
    
    res %>% 
        write_csv(here("Past Models", "Predictions", str_c("Predictions_", fitted_dates[i], ".csv")))
    
    rm(res)
    gc()
}


