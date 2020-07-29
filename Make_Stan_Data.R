Make_Stan_Data <- function(min_case_rate = 0.02, 
                           min_days = 12, 
                           min_pop = 3e5,
                           min_new_cases = 40,
                           countries = NULL, 
                           only_first_wave = TRUE,
                           stop_date = Sys.Date()) {
    options(tidyverse.quiet = TRUE)
    
    library(tidyverse, verbose = F, warn.conflicts = F, attach.required = T)
    library(readxl, verbose = F, warn.conflicts = F, attach.required = T)
    library(lubridate, verbose = F, warn.conflicts = F, attach.required = T)
    suppressPackageStartupMessages(library(mgcv, verbose = F, warn.conflicts = F, attach.required = T))
    library(broom, verbose = F, warn.conflicts = F, attach.required = T)
    
    col_types <- cols(
        .default = col_double(),
        iso_code = col_character(),
        continent = col_character(),
        location = col_character(),
        date = col_date(),
        total_tests = col_number(),
        new_tests = col_number(),
        tests_units = col_character()
    )
    
    d <- read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv",
                  col_types = col_types) %>%
        select(country = location, location = location, continent, 
               date, 
               new_cases, total_cases, 
               new_deaths, total_deaths, 
               population, population_density, median_age, gdp_per_capita, diabetes_prevalence) %>% 
        filter(!location %in% c("World", "China")) %>% 
        arrange(country, location, date) %>% 
        group_by(location) %>% 
        mutate(rolling_new_cases = data.table::frollmean(pmax(new_cases, 0), n = 7)) %>% 
        filter(any(rolling_new_cases >= min_new_cases)) %>% 
        ungroup %>% 
        mutate(case_rate = total_cases / population * 1000) %>% 
        drop_na(rolling_new_cases) %>% 
        filter(population > min_pop,
               case_rate >= min_case_rate) %>% 
        group_by(location) %>% 
        mutate(days = row_number() - 1) %>% 
        ungroup %>% 
        filter(new_cases >= 0) %>% 
        group_by(location) %>% 
        filter(n() >= min_days) %>% 
        ungroup 
    
    if (!is.null(countries)) d <- d %>% filter(country %in% countries)
    
    
    if(only_first_wave) {
        stop_time <- d %>% 
            group_by(location) %>% 
            group_nest() %>% 
            mutate(model = map(data, ~ suppressWarnings(gam(new_cases ~ s(days, bs = "gp"), 
                                                            data = ., 
                                                            family = nb(), 
                                                            method = "REML"))),
                   data = map2(data, model, 
                               function(data, model) {
                                   data$pred_new <- predict(model,type="response")
                                   data$pred_dnewdt <- c(NA, diff(data$pred_new))
                                   data$pred_log_new <- predict(model)
                                   data$pred_log_dnewdt <- c(NA, diff(data$pred_log_new))
                                   data
                               })
            ) %>% 
            select(-model) %>% 
            unnest(data) %>% 
            group_by(location) %>% 
            mutate(positive_df = 1 * (pred_log_dnewdt > 0),
                   any_negative_df = 1 * (cumsum(pred_log_dnewdt < 0 & days > 30) > 0),
                   stop = suppressWarnings(min(days[which(positive_df & any_negative_df)])),
                   stop_date = min(date) + stop) %>% 
            distinct(location, stop)
        
        d <- d %>% 
            inner_join(stop_time, by = "location") %>% 
            filter(days < stop, date <= stop_date) %>% 
            group_by(location) %>%
            filter(n() >= min_days) %>%
            ungroup %>%
            mutate(location_id = as.numeric(as.factor(location))) %>% 
            select(-stop)
        
    } else {
        d <- d %>%  
            group_by(country) %>%  
            filter(n() >= min_days, any(case_rate < upper_mult * min_case_rate)) %>% 
            ungroup %>% 
            mutate(location_id = as.numeric(as.factor(location)))
    }
    
    return(d)
}
