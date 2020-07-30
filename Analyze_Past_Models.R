library(here)
library(tidyverse)
library(broom)
library(cowplot)
library(rstan)
library(tidybayes)
library(scales)
library(lubridate)
library(plotly)

theme_set(theme_classic(base_size = 12) + 
              background_grid(minor = "none", major = "none") +
              theme(legend.position = "none"))

fitted_dates <- here("Past Models", "Models") %>% 
    list.files %>% 
    str_match("_(2020-[0-9]{2}-[0-9]{2})\\.rds") %>% 
    .[, 2] %>% 
    na.omit %>% 
    ymd

d <- here("Past Models", "Data", str_c("COVID_Data_", ymd("2020-07-25"), ".csv")) %>% 
    read_csv

make_predictions <- function(plot_locations, plot_dates) {
    
    read_fun <- function(plot_date) {
        here("Past Models", "Predictions", str_c("Predictions_", plot_date, ".csv")) %>% 
            read_csv(col_types = cols()) %>% 
            mutate(plot_date = plot_date)
    }
    
    results <- map(plot_dates, read_fun) %>% 
        reduce(bind_rows) %>% 
        filter(location %in% plot_locations,
               date <= ymd("2020-08-01"), 
               variable == "new_cases")
    
    results %>% 
        mutate(conf.high = ifelse(location == "Mexico", pmin(1e4, conf.high), conf.high)) %>% 
        ggplot(aes(date, median)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
        geom_vline(aes(xintercept = plot_date), lty = 2, alpha = 0.4) +
        geom_line() +
        geom_point(
            data = d %>%
                filter(location %in% plot_locations) %>% 
                pivot_longer(c(new_cases), names_to = "variable") %>% 
                expand_grid(plot_date = plot_dates) %>% 
                inner_join(
                    results %>%
                        group_by(location) %>% 
                        summarise(min_plot_date = min(plot_date)),
                    by = "location") %>% 
                filter(plot_date >= min_plot_date),
            aes(x = date, y = value, col = date > plot_date), size = 1
        ) +
        scale_colour_manual(values = c("black", "grey60")) +
        scale_y_continuous(breaks = pretty_breaks(3), limits = c(0, NA), expand = expansion(mult = 0.01)) +
        scale_x_date(breaks = ymd("2020-04-01", "2020-05-01", "2020-06-01", "2020-07-01"),
                     date_labels = "%B", 
                     expand = expansion(add = 10),
                     limits = c(NA_Date_, ymd("2020-08-01"))) +
        facet_grid(location ~ plot_date, scales = "free_y") +
        theme(axis.title = element_blank()) +
        background_grid(major = "none", minor = "none")
    
}



make_predictions(
    c(
        "Italy",
        "Spain",
        "Belgium",
        "Germany",
        "United Kingdom"
    ),
    ymd(c("2020-04-01", "2020-05-01", "2020-06-01"))
) +
    coord_cartesian(xlim = c(NA_Date_, ymd("2020-07-21"))) +
    ggsave(here("Past Models", "Results", "Figures", "euro_preds.png"),
           width = 6, height = 0.621 * 6, scale = 2)



make_predictions(
    c(
        "Russia",
        "Brazil",
        "Colombia",
        "Indonesia",
        "Mexico"
    ),
    ymd(c("2020-05-01", "2020-06-01", "2020-07-01"))
) +
    ggsave(here("Past Models", "Results", "Figures", "global_preds.png"),
           width = 6, height = 0.621 * 6, scale = 2)
