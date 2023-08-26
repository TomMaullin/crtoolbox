library(ggplot2)
library(readr) # for reading CSVs
library(dplyr) # for data manipulation

# Number of subjects
n <- 25

# Read in results from the appropriate directory
clt_est_store <- as.vector(read_csv(file.path(data_dir, 'clt_est_store.csv')))
sim_integral_values <- as.vector(read_csv(file.path(data_dir, 'sim_integral_values.csv')))
boot_integral_values <- read_csv(file.path(data_dir, 'boot_integral_values.csv'))

# Create data for plots
empirical_clt_df <- data.frame(Value = sort(clt_est_store), 
                               ECDF = seq(0, 1, length.out = length(clt_est_store)),
                               Label = "Empirical (muhat-mu)/tau")

empirical_integral_df <- data.frame(Value = sort(sim_integral_values), 
                                    ECDF = seq(0, 1, length.out = length(sim_integral_values)),
                                    Label = "Empirical Integral")

# Data for ribbon (shaded area)
ribbon_df <- data.frame(Value = apply(boot_integral_values, 2, median), 
                        ymin = apply(boot_integral_values, 2, quantile, probs = 0.05),
                        ymax = apply(boot_integral_values, 2, quantile, probs = 0.95),
                        ECDF = seq(0, 1, length.out = ncol(boot_integral_values)))

# Plot using ggplot2
ggplot() +
  geom_line(data = empirical_clt_df, aes(x = Value, y = ECDF, color = Label)) +
  geom_line(data = empirical_integral_df, aes(x = Value, y = ECDF, color = Label)) +
  geom_ribbon(data = ribbon_df, aes(x = Value, ymin = ymin, ymax = ymax), fill = "turquoise4", alpha = 0.5) +
  scale_color_manual(values = c("Empirical (muhat-mu)/tau" = "salmon", "Empirical Integral" = "darkorchid")) +
  labs(title = paste("Simulation results:", n, "subjects"),
       x = "Value",
       y = "ECDF") +
  theme_minimal() +
  theme(legend.position = "bottom")
