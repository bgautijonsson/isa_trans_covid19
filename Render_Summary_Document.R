library(stringr)

# Change this date to get results from other fitting dates.
fit_date <- "2020-07-28"

if (!file.exists(str_c("Past Models/Results/", fit_date))) {
    dir.create(str_c("Past Models/Results/", fit_date))
    dir.create(str_c(str_c("Past Models/Results/", fit_date), "/Figures"))
}

Sys.sleep(1)

rmarkdown::render(
    input = "Past Models/Results/Model_Parameters.Rmd",
    params = list(
        date = fit_date
    ),
    output_format = "html_document",
    output_file = str_c(fit_date, "/Model_Parameters_", fit_date, ".html")
)

