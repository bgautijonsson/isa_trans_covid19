library(readr)
library(dplyr)
library(rstan)
library(magrittr)
library(stringr)
library(lubridate)
library(here)

options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
source("Make_Stan_Data.R")

dates_to_fit <- c(ymd(c(
    "2020-04-01",
    "2020-05-01",
    "2020-06-01",
    "2020-07-01",
    "2020-07-28"
)))

countries <- c(
    "Afghanistan", "Algeria", "Argentina", "Armenia", "Austria", "Azerbaijan",
    "Bahrain", "Bangladesh", "Belarus", "Belgium", "Bolivia", "Brazil",
    "Canada", "Chile", "Colombia", "Cote d'Ivoire", "Croatia", "Czech Republic", 
    "Denmark", "Dominican Republic",
    "El Salvador", "Egypt",
    "Finland", "France", 
    "Germany", "Greece", "Guatemala",
    "Honduras", "Hungary",
    "Iceland", "India", "Indonesia", "Iran", "Iraq", "Ireland", "Israel", "Italy",
    "Japan",
    "Kenya", "Kuwait",
    "Luxembourg",
    "Mexico", 
    "Netherlands", "Nepal", "Nigeria", "Norway", 
    "Oman",
    "Pakistan", "Philippines", "Peru", "Poland", "Portugal",
    "Qatar", 
    "Russia",
    "Saudi Arabia", "Serbia", "Singapore", "Slovenia", "South Africa", 
    "Spain", "Sudan", "Sweden", "Switzerland",  
    "Turkey",   
    "Ukraine", "United Arab Emirates", "United Kingdom",
    "Venezuela"
)

for (i in seq_along(dates_to_fit)) {
    
    if (file.exists(here("Past Models", "Models", str_c("Model_", dates_to_fit[i], ".rds")))) {
        next
    }
    
    writeLines(str_c("\n",
                     "Fitting model for: ", 
                     dates_to_fit[i],
                     "\n"))
    
    
    d <- Make_Stan_Data(stop_date = dates_to_fit[i], countries = countries)
    
    N_obs <- nrow(d)
    N_locations <- max(d$location_id)
    
    
    days <- d$days
    new_cases <- d$new_cases
    location <- d$location_id %>% as.integer
    
    population <- d %>% distinct(location_id, population) %>% arrange(location_id) %>%  .$population
    stan_data <- list(N_obs = N_obs,
                      N_locations = N_locations,
                      days = days, 
                      new_cases = new_cases, 
                      location = location,
                      population = population)
    
    if (dates_to_fit[i] < ymd("2020-07-01")) {
        model_file <- "Stan/Hierarchical_InvGenLogistic_NegBin_Gompertz_MVNonCentered.stan"
    } else {
        model_file <- "Stan/Hierarchical_InvGenLogistic_NegBin_Gompertz_MVCentered.stan"
    }
    
    m <- stan(
        file = model_file, 
        data  = stan_data, 
        chains = 4, 
        iter = 2500, 
        warmup = 500,
        cores = 4,
        save_warmup = FALSE,
        control = list(max_treedepth = 15)
    )
    
    write_rds(m, here("Past Models", "Models", str_c("Model_", dates_to_fit[i], ".rds")))
    write_csv(d, here("Past Models", "Data", str_c("COVID_Data_", dates_to_fit[i], ".csv")))
}