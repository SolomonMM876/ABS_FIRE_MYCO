library(car)  # For normality tests

# Function to check normality with histograms and Shapiro-Wilk test
# Function to check normality with histograms and Shapiro-Wilk test
check_normality <- function(data, metric) {
  par(mfrow = c(2, 3))  # Arrange plots in a 2x3 grid
  
  # Compute Shapiro-Wilk tests for original and transformed data
  shapiro_original <- shapiro.test(data[[metric]])$p.value
  shapiro_log <- shapiro.test(log10(data[[metric]] + 1))$p.value
  shapiro_sqrt <- shapiro.test(sqrt(data[[metric]] + 1))$p.value
  shapiro_cuberoot <- shapiro.test((data[[metric]])^(1/3))$p.value
  shapiro_reciprocal <- shapiro.test(1 / (data[[metric]] + 1))$p.value
  shapiro_zscore <- shapiro.test(scale(data[[metric]]))$p.value
  
  # Function to add p-value text inside histogram plots
  add_p_value <- function(p_value) {
    mtext(paste("Shapiro p-value:", round(p_value, 4)), side = 3, line = -2, adj = 0.5, cex = 0.8)
  }
  
  # Generate histograms with Shapiro-Wilk p-values
  hist(data[[metric]], main = paste("Original", metric), xlab = metric)
  add_p_value(shapiro_original)
  
  hist(log10(data[[metric]] + 1), main = "Log10 Transformation", xlab = paste("log10(", metric, "+1)"))
  add_p_value(shapiro_log)
  
  hist(sqrt(data[[metric]] + 1), main = "Square Root Transformation", xlab = paste("sqrt(", metric, ")"))
  add_p_value(shapiro_sqrt)
  
  hist((data[[metric]])^(1/3), main = "Cube Root Transformation", xlab = paste("Cube Root of", metric))
  add_p_value(shapiro_cuberoot)
  
  hist(1 / (data[[metric]] + 1), main = "Reciprocal Transformation", xlab = paste("1/", metric, "+1)"))
  add_p_value(shapiro_reciprocal)
  
  hist(scale(data[[metric]]), main = "Z-score Normalization", xlab = paste("Z-score of", metric))
  add_p_value(shapiro_zscore)
  
  par(mfrow = c(1, 1))  # Reset plot layout
  
  print(paste("Shapiro-Wilk test results for", metric))
  print(data.frame(
    Transformation = c("Original", "Log10", "Square Root", "Cube Root", "Reciprocal", "Z-score"),
    P_Value = c(shapiro_original, shapiro_log, shapiro_sqrt, shapiro_cuberoot, shapiro_reciprocal, shapiro_zscore)
  ))
}


#metrics
metrics <- c('Shannon', 'Simpson', 'Chao1', 'Observed', 'Pielou')

for (metric in metrics) {
  check_normality(alpha_myco_Site_Tran, metric)#Change df here
}



