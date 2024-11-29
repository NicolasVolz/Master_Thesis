

#-------------------------------------------------------------------------------
# Simulations
#-------------------------------------------------------------------------------



## Load necessary libraries
library(parallel)
library(doParallel)
library(foreach)
library(fda)
library(MASS)
library(reshape2)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(patchwork)
library(scales)
library(pracma)
library(ggpubr)
library(grid)


#setwd("...")


################################################################################
#-----------------------DEFINE VARIABLES     -----------------------------------
################################################################################

# Sample b_values fixed for the whole simulation
set.seed(123)

# number of basis functions for the simulation
number_basis = 7    

beta_0_v <- rnorm(n = 7, mean = 3, sd = 0.1) * (1/(1:number_basis))
beta_1_v <- rnorm(n = 7, mean = 0, sd = 0)
beta_2_v <- (rnorm(n = 7, mean = -1, sd = 0.1)* (1/(1:number_basis)))/100


# Sample covariance matrices fixed for the whole simulation
Sigma_k1 <- generate_covariance_matrix (dimension = 7, df = 1.99, scale_factor = 0.1)
Sigma_k2 <- generate_covariance_matrix (dimension = 7, df = 1.99, scale_factor = 0.1)
print(round(diag(diag(Sigma_k1)),2))



### define beta  functions 
p  = length(beta_0_v)

# Define beta_0_f as a function
beta_0_f <- function(t) {
  create_beta_function(t, beta_0_v, basis)
}
# Define beta_0_f as a function
beta_1_f <- function(t, delta_2) {
  create_beta_function(t, beta_1_v, basis) + delta_2
}
# Define beta_0_f as a function
beta_2_f <- function(t) {
  create_beta_function(t, beta_2_v, basis)
}

# covariance fnctions
# arbitrary V
covariance_structure_1 = function(t_values){
  create_covariance_structure(t_values, Sigma_k1)
}

covariance_structure_2 = function(t_values){
  create_covariance_structure(t_values, Sigma_k1)
}



# diagonal V 
covariance_structure_1_1 = function(t_values){
  create_covariance_structure(t_values, diag(diag(Sigma_k1)))
}



# scaled down covariance function 
cov_2_func <- function(t_values) cov_1_func(t_values) / 5  # Scaled version for heteroscedasticity



### checks
t_values = seq(0,5,length.out=100)
check_covariance_matrix(Sigma_k1)
check_covariance_matrix(round(create_covariance_structure(t_values, Sigma_k1), digits = 6))
check_covariance_matrix(round(create_covariance_structure(t_values, diag(diag(Sigma_k1))), digits = 6))


# plot basis functions



## scale delta
Sigma_true_1 <- covariance_structure_1(t_values) #arbitrary
Sigma_true_1_1 <- covariance_structure_1_1(t_values) #diagonal

Sigma_true_2 <- covariance_structure_1(t_values)
total_variation <- (calculate_trace(Sigma_true_1) + calculate_trace(Sigma_true_2)) / 2

delta_2_print <- seq(0, 2, length.out = 10)  # Unscaled delta_2 values for clarity
delta_2_values <- delta_2_print / total_variation  # Scaled delta_2 values for simulation



## plot covariance
png(filename = "plot_V_full.png", width = 800, height = 600, res = 110)
plot_70 = plot_V_matrix(Sigma_k1, title = expression(~bold(V)))
print(plot_70)
dev.off()

png(filename = "plot_V_diagonal.png", width = 800, height = 600, res = 110)
plot_70 = plot_V_matrix(diag(diag(Sigma_k1)), title = expression(~bold(V)))
print(plot_70)
dev.off()


png(filename = "plot_covariance_diagonal.png", width = 800, height = 600, res = 110)
plot_100 = plot_covariance_matrix(Sigma_true_1_1, bquote(gamma("t", "t'")))
print(plot_100)
dev.off()

png(filename = "plot_covariance_full.png", width = 800, height = 600, res = 110)
plot_100 = plot_covariance_matrix(Sigma_true_1, bquote(gamma("t", "t'")))
print(plot_100)
dev.off()







################################################################################
#---------------------------------SIM SET UP EXAMPLE ---------------------------
################################################################################
# Number of cores to use
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

start_time <- Sys.time()

## Simulation setup
nsim = 100  # Number of simulations
results_K <- vector("numeric", length = nsim)
results_F <- vector("numeric", length = nsim)
results_L2 <- vector("numeric", length = nsim)


n1 = 30
n2 = 30
Age_1 = seq(20, 70, length.out = n1)
Age_2 = seq(20, 70, length.out = n2)

# Parallelized simulation using foreach
results <- foreach(i = 1:nsim, .combine = 'rbind', 
                   .packages = c('fda', 'MASS','reshape2','ggplot2' )) %dopar% {
                     Test = simulate_2s(n1 = n1, n2= n2,
                                        n_basis = p,
                                        beta_0_f, beta_1_f, beta_2_f,
                                        covariance_structure_1_1, covariance_structure_1_1,
                                        plot = FALSE,
                                        delta_2 = delta_2_values[1], 
                                        Age_1 = Age_1, 
                                        Age_2 = Age_2)
                     c(K = as.numeric(Test[1]), F = as.numeric(Test[2]), L2 = as.numeric(Test[3]))
                   }

# Stop the cluster
stopCluster(cl)

# Extract results
results_K <- results[, "K"]
results_F <- results[, "F"]
results_L2 <- results[, "L2"]

print(paste("power K_Test:", mean(results_K)))
print(paste("power Test F:", mean(results_F)))
print(paste("power Test L2:", mean(results_L2)))


end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(paste("Total simulation time:", elapsed_time))






  
################################################################################
#---------------------------------create Plots----------------------------------
################################################################################
#set.seed(123)
 {
  nsim = 1
  results_K <- vector("list", length = nsim)
  results_F <- vector("list", length = nsim)
  results_L2 <- vector("list", length = nsim)
  
  #delta_2 = 0.1 / (total_variation/10)
  #print(delta_2)
    
  n1 = 100
  n2 = 100
  Age_1 = seq(20, 70, length.out = n1)
  Age_2 = seq(20, 70, length.out = n2)
    
    
   # Measure the elapsed time using system.time()
    elapsed_time <- system.time({
    for (i in 1:nsim) {
      ### SIMULATION FUNCTION
      Test = simulate_2s(n1 =  n1, n2 = n2,
                           n_basis = p,
                           beta_0_f, beta_1_f, beta_2_f,
                           covariance_structure_1_1, covariance_structure_1_1,
                           plot = TRUE,
                           basis = TRUE,
                           delta_2 = delta_2_values[1], 
                           Age_1 = Age_1, 
                           Age_2 = Age_2)
    
      results_K[i] = Test[1]      
      results_F[i] = Test[2]
      results_L2[i] = Test[3]
    }
  })
  
  print(paste("power K_Test:", mean(as.numeric(results_K))))
  print(paste("power Test F:", mean(as.numeric(results_F))))
  print(paste("power Test L2:", mean(as.numeric(results_L2))))
    
  # Print elapsed time in seconds
  print(paste("Total simulation time in seconds:", elapsed_time["elapsed"]))
    
}







################################################################################
#--------------------------------- SIM -----------------------------------------
################################################################################

############################### Simulation setup ################################
library(doParallel)
library(foreach)
library(fda)
library(MASS)

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

nsim <- 20 # 10000 # Number of simulations
n_total_values <- c(20, 30, 50, 80, 120, 170, 250, 350)  # Total sample size values
delta_2_print <- seq(0, 2, length.out = 10)  # Unscaled delta_2 values for clarity
delta_2_values <- delta_2_print / total_variation  # Scaled delta_2 values for simulation
print("Delta_2 values (unscaled):")
print(delta_2_print)
print("Delta_2 values (scaled by total_variation):")
print(delta_2_values)

# Define covariance structures
covariances = list(covariance_structure_1, covariance_structure_2)

############################## Parallelized simulation ##########################
for (cov_case in c("no_covariance", "with_covariance")) {
  
  # Define the covariance functions based on the case
  if (cov_case == "no_covariance") {
    covariance_structure_1_func <- function(t_values) create_covariance_structure(t_values, diag(diag(Sigma_k1)))
  } else if (cov_case == "with_covariance") {
    covariance_structure_1_func <- function(t_values) create_covariance_structure(t_values, Sigma_k1)
  }
  
  for (c in 1:length(covariances)) {
    if (c == 1) {
      # Homoscedastic case
      homoscedastic <- TRUE
      cov_1_func <- covariance_structure_1_func
      cov_2_func <- covariance_structure_1_func  # Same function for both groups in homoscedastic case
    } else {
      # Heteroscedastic case
      homoscedastic <- FALSE
      cov_1_func <- covariance_structure_1_func
      cov_2_func <- function(t_values) cov_1_func(t_values) / 5  # Scaled version for heteroscedasticity
    }
    
    for (d in 1:length(delta_2_values)) {
      delta_2_sim <- delta_2_values[d]
      
      # Loop over conditions for n1 and n2
      for (condition in c("n2_eq_n1", "n2_eq_2n1")) {
        results_summary <- list()
        
        for (k in 1:length(n_total_values)) {
          n_total <- n_total_values[k]
          
          if (condition == "n2_eq_n1") {
            n1 <- round(n_total / 2)
            n2 <- n_total - n1  # Ensure n1 + n2 equals n_total
          } else if (condition == "n2_eq_2n1") {
            n1 <- round(n_total / 3)
            n2 <- n_total - n1  # Ensure n1 + n2 equals n_total
          }
          
          # Check if this configuration is relevant for plotting
          is_relevant <- (delta_2_sim == 0) || (n_total == 350)
          print(paste("Checking relevance for delta_2_sim =", delta_2_sim, ", n_total =", n_total, ": is_relevant =", is_relevant))  # Debugging print
          
          Age_1 = seq(20, 70, length.out = n1)
          Age_2 = seq(20, 70, length.out = n2)
          
          if (is_relevant) {
            # Parallelized simulation using foreach
            results <- foreach(i = 1:nsim, .combine = 'rbind', .packages = c('fda', 'MASS')) %dopar% {
              Test <- simulate_2s(n1 = n1, n2 = n2,
                                  n_basis = p,
                                  beta_0_f, beta_1_f, beta_2_f,
                                  cov_1_func, cov_2_func,
                                  plot = FALSE,
                                  delta_2 = delta_2_sim, 
                                  Age_1 = Age_1, 
                                  Age_2 = Age_2)
              c(K = as.numeric(Test[1]), F = as.numeric(Test[2]), L2 = as.numeric(Test[3]))
            }
            
            # Calculate power for each test
            power_K <- mean(results[, "K"])
            power_F <- mean(results[, "F"])
            power_L2 <- mean(results[, "L2"])
          } else {
            # If not relevant, set placeholders
            power_K <- 999
            power_F <- 999
            power_L2 <- 999
          }
          
          # Create a result data frame for the current configuration
          result <- data.frame(
            delta_2 = delta_2_sim,
            n1 = n1,
            n2 = n2,
            homoscedastic = homoscedastic,
            covariance_case = cov_case,
            condition = condition,
            Power_K = power_K,
            Power_F = power_F,
            Power_L2 = power_L2
          )
          
          # Append summarized results to the list
          results_summary[[k]] <- result
          
          # Print the result for the current configuration
          print(result)
        }
        
        # Combine all summarized results into one data frame for this delta_2
        results_summary_df <- do.call(rbind, results_summary)
        
        # Save the summarized results to a file using the index d
        filename <- paste0("results_summary_", cov_case, "_delta2_", d - 1, "_cov_", c, "_", condition, ".csv")
        write.csv(results_summary_df, filename, row.names = FALSE)
      }
    }
  }
}

# Stop the cluster
stopCluster(cl)








# Define aesthetics for tests
test_colors <- c("L2 Test" = "blue", "F Test" = "red", "K Test" = "orange")
test_shapes <- c("L2 Test" = 18, "F Test" = 17, "K Test" = 16)
test_line_types <- c("L2 Test" = "solid", "F Test" = "dashed", "K Test" = "dotted")

scenario_titles <- c(
  "Bal. Homo. Uncorr.",
  "Bal. Hetero. Uncorr.",
  "Unbal. Hetero. Uncorr.",
  "Unbal. Homo. Uncorr.",
  "Bal. Homo. Corr.",
  "Bal. Hetero. Corr.",
  "Unbal. Hetero. Corr.",
  "Unbal. Homo. Corr."
)



load_data <- function(delta_indices, file_prefix, sample_filter = NULL) {
  combined_data <- data.frame()
  
  for (cov_case in c("no_covariance", "with_covariance")) {
    for (c in 1:2) {
      homoscedastic <- ifelse(c == 1, "Homo.", "Hetero.")
      for (condition in c("n2_eq_n1", "n2_eq_2n1")) {
        balance <- ifelse(condition == "n2_eq_n1", "Bal.", "Unbal.")
        correlation <- ifelse(cov_case == "no_covariance", "Uncorr.", "Corr.")
        
        # Set the scenario title based on attributes in the desired order
        scenario_title <- paste(balance, homoscedastic, correlation)
        
        for (d in delta_indices) {
          # Use the index directly in the filename
          filename <- paste0(file_prefix, cov_case, "_delta2_", d, "_cov_", c, "_", condition, ".csv")
          
          if (file.exists(filename)) {
            results_summary_df <- read.csv(filename)
            
            # Apply sample size filter if provided
            if (!is.null(sample_filter)) {
              results_summary_df <- subset(results_summary_df, n1 + n2 == sample_filter)
            }
            
            # Melt the data for ggplot2
            df_melted <- melt(results_summary_df, id.vars = c("delta_2", "n1", "n2", "homoscedastic", "covariance_case", "condition"), 
                              measure.vars = c("Power_L2", "Power_F", "Power_K"), variable.name = "Test", value.name = "Power")
            
            # Rename test types for readability
            df_melted$Test <- factor(df_melted$Test, levels = c("Power_L2", "Power_F", "Power_K"), labels = c("L2 Test", "F Test", "K Test"))
            
            # Add dynamically generated scenario title and a new "Scenario_Order" for ordering by the desired layout
            df_melted$Scenario <- factor(scenario_title, levels = scenario_titles)
            
            # Combine with the main data frame
            combined_data <- rbind(combined_data, df_melted)
          } else {
            print(paste("File not found:", filename))
          }
        }
      }
    }
  }
  return(combined_data)
}

new_titles <- c(
  "Setting 1",
  "Setting 2",
  "Setting 3",
  "Setting 4",
  "Setting 5",
  "Setting 6",
  "Setting 7",
  "Setting 8"
)

# Load data for Type 1 Error plot using indices
type1_data <- load_data(delta_indices = c(0), file_prefix = "results_summary_", sample_filter = NULL)

# Reassign the Scenario factor with proper levels and labels
type1_data$Scenario <- factor(
  type1_data$Scenario,
  levels = scenario_titles,  # Original scenario order
  labels = new_titles       # Desired titles
)

# Define the indices for delta_2_values as a sequence from 0 to the length of delta_2_values - 1
delta_indices <- 0:(length(delta_2_values) - 1)

# Load the data using the indices as delta_values
power_data <- load_data(delta_indices = delta_indices, file_prefix = "results_summary_", sample_filter = 350)

# Reassign the Scenario factor with proper levels and labels
power_data$Scenario <- factor(
  power_data$Scenario,
  levels = scenario_titles,  # Original scenario order
  labels = new_titles       # Desired titles
)

# Define the plotting function with correct labels and manual facet ordering by Scenario
plot_simulation <- function(data, plot_title, x_label, y_label, x_variable = "delta_2", add_hline = FALSE, x_axis_format = NULL) {
  p <- ggplot(data, aes_string(x = x_variable, y = "Power", color = "Test", shape = "Test")) +
    geom_line(aes(linetype = Test), size = 1) +
    geom_point(size = 3) +
    labs(title = plot_title, x = x_label, y = y_label) +
    scale_color_manual(values = test_colors, name = NULL) +
    scale_shape_manual(values = test_shapes, name = NULL) +
    scale_linetype_manual(values = test_line_types, name = NULL) +
    facet_wrap(~ Scenario, ncol = 4) +
    ylim(0, 1) +  # Set y-axis limits from 0 to 1
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 18),
      strip.text = element_text(size = 16)  # Adjust the size of the facet titles here
    )
  
  # Apply x-axis formatting if specified
  if (!is.null(x_axis_format)) {
    p <- p + scale_x_continuous(labels = x_axis_format)
  }
  
  # Add a dashed horizontal line at y = 0.05 if specified
  if (add_hline) {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dashed", color = "black")
  }
  
  return(p)
}
# Generate Type 1 Error Plot with horizontal line at y = 0.05
type1_error_plot <- plot_simulation(
  data = type1_data,
  plot_title = "Type 1 Error Simulation",
  x_label = "Total Sample Size",
  y_label = "Type 1 Error",
  x_variable = "n1 + n2",
  add_hline = TRUE,
  x_axis_format = number_format(accuracy = 1)  # Integer format for Type 1 Error plot
)

# Save Type 1 Error Plot
ggsave("Combined_Type1Error_Simulation_Plot.jpg", plot = type1_error_plot, width = 14, height = 9, units = "in", dpi = 300)

# Generate Power Simulation Plot
power_simulation_plot <- plot_simulation(
  data = power_data,
  plot_title = "Power Simulation",
  x_label = expression(delta[scaled]),
  y_label = "Power",
  x_axis_format = number_format(accuracy = 0.01)  # Two decimal places for Power plot
)

# Save Power Simulation Plot
ggsave("Combined_Power_Simulation_Plot.jpg", plot = power_simulation_plot, width = 14, height = 9, units = "in", dpi = 300)





################################################################################
### create basis plots

# Fourier Basis Functions
# Create a Fourier basis
fourier_basis <- create.fourier.basis(rangeval = c(0, 2*pi), nbasis = 10)
# Evaluate the basis functions at a grid of points
x_fourier <- seq(0, 2*pi, length.out = 100)
fourier_values <- eval.basis(x_fourier, fourier_basis)




#setwd("...")
# Generate the plot (this will automatically be saved to the file)
png(filename = "Basis_Plot.png", width = 800, height = 600, res = 110)

# Create plots for each Fourier basis function with labels
plots_fourier <- lapply(1:6, function(i) plot_single_basis(x_fourier, fourier_values[, i], bquote(B[.(i-1)](t))))

# Arrange Fourier basis plots in a grid
grid.arrange(grobs = plots_fourier, ncol = 3, nrow = 2, top = textGrob("Fourier Basis Functions", gp = gpar(fontsize = 16)))

# Close the device
dev.off()




