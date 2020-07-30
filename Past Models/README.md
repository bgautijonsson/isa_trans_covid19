* **Data:** The data used to train the models is stored here after model training.
* **Models:** Each model is saved as a binary `.Rds` object and read into memory for prediction and processing when needed
* **Predictions:** The file `Predict_Past_Models.R` in the top level directory makes specific model predictions and stores them here as `.csv` files that are more easily handled than the full `.Rds` model binaries.
* **Results:** Here are files for creating an `.html` report for estimated model parameters and also a Shiny app for interactive model checking and prediction.