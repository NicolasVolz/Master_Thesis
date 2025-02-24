
# Libraries
library(tidyverse)
library(readxl)
library(parallel)
library(doParallel)
library(foreach)
library(fda)
library(MASS)
library(reshape2)
library(ggplot2)
library(dplyr)
library(corrplot)
library(gridExtra)
library(patchwork)


################################################################################
######----main regression function--------######################################
################################################################################

do_regression <- function(pain_data, nopain_data, n_pain, n_nopain, id, to_merge, evaluation_grid, model_name, kroff = FALSE, plot = FALSE) {
  
  # Step 1: Build X data frame and filter out missing age/gender
  X <- rep(1, (n_pain + n_nopain))
  pain <- c(rep(1, n_pain), rep(0, n_nopain))
  X_null <- data.frame(id, X, pain)
  X_complete <- left_join(X_null, to_merge, by = "id")
  missing_indices <- which(is.na(X_complete$age) | is.na(X_complete$gender) | X_complete$gender == 7)
  X_complete <- X_complete[-missing_indices, ]
  X_complete$gender <- ifelse(X_complete$gender == 1, 0, 1)
  n <- n_pain + n_nopain
  
  # Step 1a: Compute summary statistics for "Pain" and "No Pain"
  pain_group <- X_complete[X_complete$pain == 1, ]
  no_pain_group <- X_complete[X_complete$pain == 0, ]
  
  summary_statistics <- list(
    Pain = list(
      Gender = table(pain_group$gender),
      Age = list(
        Mean = mean(pain_group$age, na.rm = TRUE),
        SD = sd(pain_group$age, na.rm = TRUE),
        Range = range(pain_group$age, na.rm = TRUE)
      )
    ),
    No_Pain = list(
      Gender = table(no_pain_group$gender),
      Age = list(
        Mean = mean(no_pain_group$age, na.rm = TRUE),
        SD = sd(no_pain_group$age, na.rm = TRUE),
        Range = range(no_pain_group$age, na.rm = TRUE)
      )
    )
  )
  
  # Print or log the summary statistics
  print("Summary Statistics:")
  print(summary_statistics)
  
  # Step 2: Combine coefficients and create functional data objects
  data_list <- list(pain_data, nopain_data)
  combined_coefs_matrix <- combine_coefs_into_one(data_list)
  rangeval <- c(0, 1440)
  nbasis <- 25
  basisobj <- create.fourier.basis(rangeval, nbasis)
  fd_obj <- fd(coef = combined_coefs_matrix, basisobj = basisobj)
  fd_obj_complete <- fd(coef = combined_coefs_matrix[, -missing_indices], basisobj = basisobj)
  
  # Step 3: Set up indices and groups for plotting
  coef_matrix <- fd_obj_complete$coefs
  basis <- fd_obj_complete$basis
  pain_indicator <- X_complete$pain
  gender_indicator <- X_complete$gender
  pain_indices <- which(pain_indicator == 1)
  nopain_indices <- which(pain_indicator == 0)
  male_indices <- which(gender_indicator == 0)
  female_indices <- which(gender_indicator == 1)
  
  # Step 4: Create fd objects for each group
  fd <- fd(coef = coef_matrix, basis = basis)
  fd_pain <- fd(coef = coef_matrix[, pain_indices], basis = basis)
  fd_nopain <- fd(coef = coef_matrix[, nopain_indices], basis = basis)
  fd_male <- fd(coef = coef_matrix[, male_indices], basis = basis)
  fd_female <- fd(coef = coef_matrix[, female_indices], basis = basis)
  mean_fd_pain <- mean.fd(fd_pain)
  mean_fd_nopain <- mean.fd(fd_nopain)
  mean_fd_male <- mean.fd(fd_male)
  mean_fd_female <- mean.fd(fd_female)
  
  # Step 5: Plot the data
  plot_pain_vs_no_pain(fd_pain, fd_nopain, mean_fd_pain, mean_fd_nopain, "Pain vs. No Pain")
  plot_pain_vs_no_pain(fd_male, fd_female, mean_fd_male, mean_fd_female, "Gender Differences")
  
  # Step 6: Set up regression with beta list and group list
  betafdPar <- fdPar(basisobj)
  p_interaction <- 4
  betaList <- vector("list", p_interaction)
  for (j in 1:p_interaction) {
    betaList[[j]] <- betafdPar
  }
  groupList <- vector("list", p_interaction)
  groupList[[1]] <- setNames(X_complete$X, colnames(fd_obj_complete$coefs))
  groupList[[2]] <- if (!kroff) setNames(X_complete$pain, colnames(fd_obj_complete$coefs)) else setNames(X_complete$kroff, colnames(fd_obj_complete$coefs))
  groupList[[3]] <- scale(setNames(X_complete$age, colnames(fd_obj_complete$coefs)), center = TRUE, scale = FALSE)
  groupList[[4]] <- setNames(X_complete$gender, colnames(fd_obj_complete$coefs))
  names(groupList) <- c('Intercept', 'pain', 'age', 'gender')
  
  
  # Convert age and gender to numeric
  groupList[[3]] <- as.numeric(groupList[[3]])  # Age
  groupList[[4]] <- as.numeric(groupList[[4]])  # Gender
  
  
  # Step 7: Run functional regression
  fRegressList <- fRegress(fd_obj_complete, groupList, betaList)
  betaestlist <- fRegressList$betaestlist
  fitted_values <- fRegressList$yhatfdobj
  residuals_hat <- fd_obj_complete - fitted_values
  beta_intercept <- betaestlist[[1]]$fd
  beta_pain <- betaestlist[[2]]$fd
  beta_age <- betaestlist[[3]]$fd
  beta_gender <- betaestlist[[4]]$fd
  
  
  
  
  # Step 8: Compute confidence intervals and plots for residuals
  residual_values <- t(as.matrix(residuals_hat$coefs))
  gamma_hat <- var.fd(residuals_hat)
  X_matrix <- as.matrix(X_complete[, -1])
  
  if(kroff){
    ## rearrange to kroff 
    X_matrix[,2] <- X_matrix[,5] 
    X_matrix <- X_matrix[,-5]
  }
  
  
  
  if(plot){
    ## check covariance structure 
    # Generate the plot (this will automatically be saved to the file)
    gamma_hat_matrix <- cov(residual_values)
    png(filename = paste0("Covariance_", model_name, ".png"), width = 800, height = 600, res = 100)
    plot_V_matrix(gamma_hat_matrix, title = expression(~bold(hat(V))))
    # Close the device
    dev.off()
    
    ## check covariance structure 
    # Generate the plot (this will automatically be saved to the file)
    covariance <- create_covariance_structure(evaluation_grid, gamma_hat_matrix)
    png(filename = paste0("correlation_function_", model_name, ".png"), width = 800, height = 600, res = 100)
    plot_covariance_matrix(cov2cor(covariance), title = bquote(cor(t, t * "'")))
    # Close the device
    dev.off()
    
    
    ## check homoscedasticity
    residual_values_1_pain_0= residual_values[groupList$pain==0,]
    residual_values_1_pain_1 = residual_values[groupList$pain==1,]
    
    gamma_hat_matrix_1_pain_0 <- cov(residual_values_1_pain_0)
    gamma_hat_matrix_1_pain_1 <- cov(residual_values_1_pain_1)
    
    # Generate the plot (this will automatically be saved to the file)
    png(filename = paste0("V_matrix_no_pain", model_name, ".png"), width = 800, height = 600, res = 100)
    plot_V_matrix(gamma_hat_matrix_1_pain_0, title = expression(~bold(hat(V))))
    # Close the device
    dev.off()
    
    
    # Generate the plot (this will automatically be saved to the file)
    png(filename = paste0("V_matrix_pain", model_name, ".png"), width = 800, height = 600, res = 100)
    plot_V_matrix(gamma_hat_matrix_1_pain_1, title = expression(~bold(hat(V))))
    # Close the device
    dev.off()
  }
  
  
  
  
  
  
  
  
  
  XtX <- t(X_matrix) %*% X_matrix
  XtX_inv <- solve(XtX)
  XtX_inv_diag <- diag(XtX_inv)
  se_beta_intercept <- pointwise_sd(beta_intercept, XtX_inv_diag[1], gamma_hat, evaluation_grid)
  se_beta_pain <- pointwise_sd(beta_pain, XtX_inv_diag[2], gamma_hat, evaluation_grid)
  se_beta_age <- pointwise_sd(beta_age, XtX_inv_diag[3], gamma_hat, evaluation_grid)
  se_beta_gender <- pointwise_sd(beta_gender, XtX_inv_diag[4], gamma_hat, evaluation_grid)
  
  # Step 9: Plot beta functions with confidence intervals
  plot_with_ci(beta_intercept, se_beta_intercept, "Intercept", "blue", evaluation_grid)
  plot_with_ci(beta_pain, se_beta_pain, "Pain", "blue", evaluation_grid)
  plot_with_ci_age(beta_age, se_beta_age, "Age", "blue", evaluation_grid)
  plot_with_ci(beta_gender, se_beta_gender, "Gender", "blue", evaluation_grid)
  
  
  export_plot_with_ci(paste0("Beta_Function_for_Intercept_", model_name, ".png"), beta_intercept, se_beta_intercept, "Intercept", "blue", evaluation_grid)
  export_plot_with_ci(paste0("Beta_Function_for_Pain_", model_name, ".png"), beta_pain, se_beta_pain, "Pain", "blue", evaluation_grid)
  export_plot_with_ci_age(paste0("Beta_Function_for_Age_", model_name, ".png"), beta_age, se_beta_age, "Age", "blue", evaluation_grid)
  export_plot_with_ci(paste0("Beta_Function_for_Gender_", model_name, ".png"), beta_gender, se_beta_gender, "Gender", "blue", evaluation_grid)
  
  
  
  # Set up layout for individual combined plots
  par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.2)
  
  # Export Plot 1: Pain vs. No Pain for Work Days
  png(paste0("Pain_vs_No_Pain_", model_name, ".png"), width = 1000, height = 800, res = 120)
  par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.2)
  legend_labels <- c("Pain (Mean)", "No Pain (Mean)", "Pain (Data)", "No Pain (Data)")
  plot_pain_vs_no_pain(fd_pain, fd_nopain, mean_fd_pain, mean_fd_nopain, "Pain vs. No Pain", legend_labels)
  plot_with_ci(beta_pain, se_beta_pain, "", "blue", evaluation_grid)
  dev.off()
  
  # Export Plot 2: Gender Differences for Work Days
  png(paste0("Gender_Differences_", model_name, ".png"), width = 1000, height = 800, res = 120)
  par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.2)
  legend_labels <- c("Male (Mean)", "Female (Mean)", "Male (Data)", "Female (Data)")
  plot_pain_vs_no_pain(fd_male, fd_female, mean_fd_male, mean_fd_female, "Gender Differences", legend_labels)
  plot_with_ci(beta_gender, se_beta_gender, "", "blue", evaluation_grid)
  dev.off()
  
  
  # Reset layout to default
  par(mfrow = c(1, 1))
  
  
  
  
  
  # Step 10: Additional regression testing
  n_basis <- 25
  pain_hat <- c(beta_pain$coefs)
  age_hat <- c(beta_age$coefs)
  gender_hat <- c(beta_gender$coefs)
  results <- perform_tests(
    beta_list = list(pain_hat, age_hat, gender_hat), 
    residual_values = residual_values, 
    XtX_inv_diag = XtX_inv_diag, 
    n_basis = n_basis, 
    n = n, 
    p_interaction = 3
  )
  
  # Return results for further processing
  return(list(
    fd, regression_results = results
  ))
}





do_intercept_only_regression <- function(coef_matrix_dif, evaluation_grid, model_name) {
  # Step 1: Create basis and fd object for the intercept-only model
  rangeval <- c(0, 1440)
  nbasis <- 25
  basisobj <- create.fourier.basis(rangeval, nbasis)
  fd_obj_dif <- fd(coef = coef_matrix_dif, basisobj = basisobj)
  
  # Step 2: Set up regression for intercept-only model
  X <- rep(1, ncol(coef_matrix_dif))
  groupList_intercept <- list(Intercept = setNames(X, colnames(fd_obj_dif$coefs)))
  betafdPar <- fdPar(basisobj)
  betaList_intercept <- list(betafdPar)
  
  # Run functional regression
  fRegressList_intercept <- fRegress(fd_obj_dif, groupList_intercept, betaList_intercept)
  beta_intercept_only <- fRegressList_intercept$betaestlist[[1]]$fd
  residuals_hat_intercept <- fd_obj_dif - fRegressList_intercept$yhatfdobj
  
  # Step 3: Plot the overall difference functions
  plot(fd_obj_dif, main = paste("Functional Data with Intercept Only -", model_name))
  lines(beta_intercept_only, lwd = 3, col = "black")
  
  # Step 4: Residuals and covariance for intercept-only model
  residual_values_intercept <- t(as.matrix(residuals_hat_intercept$coefs))
  gamma_hat_intercept <- var.fd(residuals_hat_intercept)
  XtX_intercept <- sum(X^2)
  XtX_inv_diag_intercept <- 1 / XtX_intercept
  
  # Step 5: Compute confidence intervals for the intercept
  se_beta_intercept_only <- pointwise_sd(beta_intercept_only, XtX_inv_diag_intercept, gamma_hat_intercept, evaluation_grid)
  
  
  # Export Plot 3: Pain vs. No Pain for Non-Work Days
  png("Work_vs_no_work.png", width = 1000, height = 800, res = 120)
  par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.2)
  legend_labels <- c("Work (Mean)", "No work (Mean)", "Work (Data)", "No work (Data)")
  plot_pain_vs_no_pain(fd_work, fd_nonwork, mean_fd_work, mean_fd_nonwork, "Work vs. No work", legend_labels)
  plot_with_ci(beta_intercept_only, se_beta_intercept_only, "", "blue", evaluation_grid)
  dev.off()
  
  
  '  # Export Plot 1: Pain vs. No Pain for Work Days
  png(paste0("Work_vs_no_work.png"), width = 1000, height = 800, res = 120)
  par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.2)
  legend_labels <- c("Work (Mean)", "No work (Mean)", "Work (Data)", "No work (Data)")
  plot_pain_vs_no_pain(fd_work, fd_nonwork, mean_fd_work, mean_fd_nonwork, "Work vs. No work", legend_labels)
  plot_with_ci(beta_intercept_only, se_beta_intercept_only, "", "blue", evaluation_grid)
  dev.off()'
  
  
  
  # Step 6: Perform K, L2, and F tests for intercept
  intercept_only_hat <- c(beta_intercept_only$coefs)
  V_hat_intercept <- cov(residual_values_intercept)
  K_intercept <- c(c(intercept_only_hat) %*% solve(V_hat_intercept) / XtX_inv_diag_intercept) %*% c(intercept_only_hat)
  significant_K_intercept <- K_intercept > qchisq(0.95, df = nbasis)
  p_value_intercept <- 1 - pchisq(K_intercept, df = nbasis)
  
  # L2 Test
  Cov_intercept <- V_hat_intercept * XtX_inv_diag_intercept
  svd_intercept <- svd(Cov_intercept)
  eigenvalues_intercept <- svd_intercept$d
  chi_squares_intercept <- sapply(eigenvalues_intercept, function(lambda) {
    rchisq(1, df = 1) * lambda
  })
  
  gamma_otimes2_intercept <- compute_gamma_otimes2(Cov_intercept)
  trace_gamma_intercept <- sum(diag(Cov_intercept))
  trace_gamma_otimes2_intercept <- sum(diag(gamma_otimes2_intercept))
  
  beta_hat_k_intercept <- trace_gamma_otimes2_intercept / trace_gamma_intercept
  k_hat_k_intercept <- (trace_gamma_intercept^2) / trace_gamma_otimes2_intercept
  critical_value_2_intercept <- beta_hat_k_intercept * qchisq(0.95, df = k_hat_k_intercept)
  
  Tn_K_actual_intercept <- (t(intercept_only_hat) %*% intercept_only_hat)
  significant_L2_intercept <- Tn_K_actual_intercept > critical_value_2_intercept
  p_value_L2_intercept <- 1 - pchisq(Tn_K_actual_intercept / beta_hat_k_intercept, df = k_hat_k_intercept)
  
  # F Test
  Fn_intercept <- Tn_K_actual_intercept / trace_gamma_intercept
  df1_k_intercept <- k_hat_k_intercept
  df2_k_intercept <- (ncol(coef_matrix_dif) - 1) * k_hat_k_intercept
  critical_value_3_intercept <- qf(0.95, df1 = df1_k_intercept, df2 = df2_k_intercept)
  significant_F_k_intercept <- Fn_intercept > critical_value_3_intercept
  p_value_F_intercept <- 1 - pf(Fn_intercept, df1 = df1_k_intercept, df2 = df2_k_intercept)
  
  # Return test results for further processing
  return(list(
    K_test = list(value = K_intercept, significant = significant_K_intercept, p_value = p_value_intercept),
    L2_test = list(value = Tn_K_actual_intercept, significant = significant_L2_intercept, p_value = p_value_L2_intercept),
    F_test = list(value = Fn_intercept, significant = significant_F_k_intercept, p_value = p_value_F_intercept)
  ))
}





################################################################################
###########--- Helper Functions----#############################################
################################################################################


###  Read Data
read_and_transform_data <- function(file_path, sep = ";", dec = ",", xlsx = TRUE, variable = "MET") {
  if(xlsx == FALSE){
  data_raw <- read.csv(file_path, sep = sep, header = TRUE, dec = dec)}
  
  if(xlsx == TRUE){
    data_raw <- read_excel(file_path)}
    
  column_names <- c("Time_rel_s", "Day_rel_d", "Time_rel_hh_mm_ss", "Date_abs_yyyy_mm_dd", "Time_abs_hh_mm_ss",
                    "Marker_Time_Type_Comment", "ActivityClass", "ActivityEnergyExpenditure_kcal_d", "Altitude_m",
                    "BodyPosition", "InclinationDown_deg", "InclinationForward_deg", "InclinationRight_deg", "MET",
                    "MovementAcceleration_g", "NonWearTime", "StepCount_steps", "TempMean", "TotalEnergyExpenditure_kcal_d",
                    "VerticalSpeed_m_s", "charging", "movementAcceleration_g", "stateOfCharge_percent", "stepCount_steps",
                    "tempMean")
  colnames(data_raw) <- column_names
  
  data_filtered <- data_raw %>%
    filter(Day_rel_d < 13)
  
  data <- data.frame(
    t_values = rep(1:1440, times = 13),
    Day_rel_d = rep(1:13, each = 1440),
    MovementAcceleration_g = data_filtered[[variable]]
  )
  
  data_wide <- data %>%
    pivot_wider(
      names_from = Day_rel_d,
      values_from = MovementAcceleration_g,
      names_prefix = "MovementAcceleration_g_d"
    ) %>%
    mutate(across(starts_with("MovementAcceleration_g_d"), as.numeric)) %>%
    na.omit()
  
  return(data_wide)
}




###  Normalize Data
normalize_data <- function(data, threshold = 0.01) {
  normalized_data <- data
  
  for (day in 2:ncol(data)) {
    values <- data[[day]]
    
    # Detect long intervals of zero values
    zero_intervals <- rle(values < threshold)
    zero_lengths <- zero_intervals$lengths
    zero_values <- zero_intervals$values
    
    # Find the longest zero interval
    if (any(zero_values & zero_lengths >= 20)) {
      start_idx <- which(zero_values & zero_lengths >= 180)[1]
      start_pos <- sum(zero_lengths[1:(start_idx - 1)]) + 1
      end_pos <- start_pos + zero_lengths[start_idx] - 1
      
      # Remove the zero interval and concatenate remaining values
      non_zero_values <- c(values[(end_pos + 1):length(values)], values[1:(start_pos - 1)])
      non_zero_times <- seq_along(non_zero_values)
      
      # Interpolate to create 1440 points from non-zero interval
      interpolated_values <- approx(non_zero_times, non_zero_values, n = 1440)$y
      
      normalized_data[[day]] <- interpolated_values
    } else {
      normalized_data[[day]] <- approx(seq_along(values), values, n = 1440)$y
    }
  }
  
  return(normalized_data)
}
    



###  Smooth Data
smooth_data <- function(sim_data, nbasis = 6, basis_type = "fourier") {
  if (is.data.frame(sim_data)) {
    day <- sim_data[[1]]
    samples_matrix <- as.matrix(sim_data[, -1])
  } else {
    day <- seq_along(sim_data)
    samples_matrix <- matrix(sim_data, ncol = 1)
  }
  
  basis <- switch(basis_type,
                  "fourier" = create.fourier.basis(c(0, max(day)), nbasis),
                  "bspline" = create.bspline.basis(c(0, max(day)), nbasis),
                  stop("Invalid basis_type. Please choose 'fourier' or 'bspline'."))
  
  smooth_result <- smooth.basis(day, samples_matrix, basis)
  smoothed_fd <- smooth_result$fd
  
  return(list(fd = smoothed_fd, gcv = smooth_result$gcv))
}





## check function 
# Function to read a single file into a data frame
read_data_file <- function(file_path) {
  tryCatch({
    data <- read_excel(file_path)
    return(data)
  }, error = function(e) {
    message(paste("Error reading file:", file_path, "\n", e))
    return(NULL)
  })
}


create_overall_mean_function <- function(results) {
  # Extract mean functions and combine them into a single data frame
  all_mean_functions <- do.call(cbind, lapply(results, function(res) res$mean_function_raw))
  
  # Calculate the overall mean function
  overall_mean_function <- rowMeans(all_mean_functions, na.rm = TRUE)
  
  # Smooth the overall mean function
  smoothed_overall_mean_fd <- smooth_data(overall_mean_function, nbasis = 25, basis_type = "fourier")
  
  return(list(overall_mean_function = overall_mean_function, smoothed_overall_mean_fd = smoothed_overall_mean_fd))
}








# Function to filter out results with five consecutive identical measurements (up to 3 decimal places)
# but only if those measurements are greater than 0.1
filter_results <- function(results) {
  # Create an empty list to store filtered results
  filtered_results <- list()
  
  # Iterate over each ID in the results list
  for (id in names(results)) {
    # Extract the mean_function_raw for the current ID
    mean_function_raw <- results[[id]]$mean_function_raw
    
    # Check if mean_function_raw is not NULL and is a numeric vector
    if (is.null(mean_function_raw) || !is.numeric(mean_function_raw)) {
      next  # Skip to the next ID if the mean_function_raw is not valid
    }
    
    # Round the mean function values to 3 decimal places
    rounded_values <- round(mean_function_raw, 5)
    
    # Find runs of identical values
    runs <- rle(rounded_values)
    
    # Check if there are at least 5 consecutive identical values greater than 0.1
    identical_sequence_found <- any(runs$lengths >= 10 & runs$values > 0.1)
    
    # If no such sequence is found, add the result to filtered_results
    if (!identical_sequence_found) {
      filtered_results[[id]] <- results[[id]]
    }
  }
  
  return(filtered_results)
}



# Function to filter out results with five consecutive identical measurements (up to 3 decimal places)
# and record which IDs were excluded with their violation time points
filter_results_with_violations <- function(results) {
  # Create a list to store filtered results and a list to track exclusions with violation times
  filtered_results <- list()
  exclusions <- list()  # To store IDs with specific time points causing exclusion
  
  # Iterate over each ID in the results list
  for (id in names(results)) {
    # Extract the mean_function_raw for the current ID
    mean_function_raw <- results[[id]]$mean_function_raw
    
    # Check if mean_function_raw is not NULL and is a numeric vector
    if (is.null(mean_function_raw) || !is.numeric(mean_function_raw)) {
      next  # Skip to the next ID if the mean_function_raw is not valid
    }
    
    # Round the mean function values to 3 decimal places
    rounded_values <- round(mean_function_raw, 5)
    
    # Find runs of identical values
    runs <- rle(rounded_values)
    
    # Identify any sequences of 10 or more consecutive identical values > 0.1
    violating_sequences <- which(runs$lengths >= 10 & runs$values > 0.1)
    
    if (length(violating_sequences) > 0) {
      # Track the start time points for each violating sequence
      violation_times <- c()  # Store indices of violations
      cumulative_index <- 0  # Track cumulative index for each run
      
      # Calculate start and end indices of each violating sequence
      for (seq_index in violating_sequences) {
        start_index <- cumulative_index + 1
        end_index <- start_index + runs$lengths[seq_index] - 1
        violation_times <- c(violation_times, start_index:end_index)
        
        # Update cumulative index
        cumulative_index <- cumulative_index + runs$lengths[seq_index]
      }
      
      # Store the violation times for this excluded ID
      exclusions[[id]] <- violation_times
    } else {
      # If no violation is found, add the result to filtered_results
      filtered_results[[id]] <- results[[id]]
    }
  }
  
  return(list(filtered_results = filtered_results, exclusions = exclusions))
}




# Print the IDs with exclusion reasons in a more readable way
print_exclusions <- function(exclusions, category) {
  cat("Excluded IDs with violations for", category, ":\n")
  if (length(exclusions) == 0) {
    cat("No exclusions found.\n")
  } else {
    for (id in names(exclusions)) {
      cat("ID:", id, "Violations at time points:", exclusions[[id]], "\n")
    }
  }
}




# Extract the IDs from the file names
extract_id <- function(filename) {
  # Assuming the ID is the part of the filename after the first underscore and before "_Results.xlsx"
  id <- sub(".*_(u[0-9]+)_Results\\.xlsx", "\\1", filename)
  return(id)
}


# Function to standardize ID format with leading zeros for numbers under 1000
standardize_id <- function(id) {
  prefix <- sub("[0-9]+", "", id)  # Extract the prefix (non-numeric part)
  numeric_part <- as.numeric(sub("[^0-9]", "", id))  # Extract the numeric part and convert to numeric
  formatted_id <- paste0(prefix, sprintf("%03d", numeric_part))  # Add leading zeros
  return(formatted_id)
}


# Function to format IDs with leading zeros for numbers under 100
format_id <- function(id) {
  prefix <- sub("[0-9]+", "", id)  # Extract the prefix (non-numeric part)
  numeric_part <- as.numeric(sub("[^0-9]", "", id))  # Extract the numeric part and convert to numeric
  if (numeric_part < 100) {
    formatted_id <- paste0(prefix, sprintf("%03d", numeric_part))  # Add leading zeros
  } else {
    formatted_id <- id  # Keep as is if 100 or above
  }
  return(formatted_id)
}




# Function to combine smoothed mean functions into a data frame with a label
combine_smoothed_mean_functions <- function(smoothed_mean_fd, label) {
  data.frame(
    t_values = 1:1440,
    mean_function = eval.fd(1:1440, smoothed_mean_fd),
    label = label
  )
}




# Function to split data based on work status
split_data_by_work_status <- function(data, work_status) {
  work_days <- which(work_status == 1)
  non_work_days <- which(work_status == 0)
  
  work_columns <- c("t_values", paste0("MovementAcceleration_g_d", work_days))
  non_work_columns <- c("t_values", paste0("MovementAcceleration_g_d", non_work_days))
  
  work_data <- data[, work_columns, drop = FALSE]
  non_work_data <- data[, non_work_columns, drop = FALSE]
  
  list(work_data = work_data, non_work_data = non_work_data)
}


# Function to process each file and split into work/non-work
process_file_1 <- function(file_path, work_status) {
  # Read and transform data
  data_wide <- read_and_transform_data(file_path)
  # Split data into work and non-work days
  split_data <- split_data_by_work_status_1(data_wide, work_status)
  work_data <- split_data$work_data
  non_work_data <- split_data$non_work_data
  
  results <- list(work = NULL, non_work = NULL)
  
  if (ncol(work_data) > 1) { # Ensure there are work day columns to process
    # Smooth data for work days
    smoothed_data_work <- smooth_data(work_data, nbasis = 25, basis_type = "fourier")
    
    # Calculate mean function for work days
    mean_function_work <- calculate_mean_function(work_data)
    mean_function_df_work <- data.frame(t_values = work_data$t_values, mean_function = mean_function_work)
    
    # Smooth mean function for work days
    smoothed_mean_fd_work <- smooth_data(mean_function_work, nbasis = 25, basis_type = "fourier")
    
    results$work <- list(
      mean_function_raw = mean_function_work,
      smoothed_mean_fd = smoothed_mean_fd_work,
      mean_function_df = mean_function_df_work
    )
  }
  
  if (ncol(non_work_data) > 1) { # Ensure there are non-work day columns to process
    # Smooth data for non-work days
    smoothed_data_non_work <- smooth_data(non_work_data, nbasis = 25, basis_type = "fourier")
    
    # Calculate mean function for non-work days
    mean_function_non_work <- calculate_mean_function(non_work_data)
    mean_function_df_non_work <- data.frame(t_values = non_work_data$t_values, mean_function = mean_function_non_work)
    
    # Smooth mean function for non-work days
    smoothed_mean_fd_non_work <- smooth_data(mean_function_non_work, nbasis = 25, basis_type = "fourier")
    
    results$non_work <- list(
      mean_function_raw = mean_function_non_work,
      smoothed_mean_fd = smoothed_mean_fd_non_work,
      mean_function_df = mean_function_df_non_work
    )
  }
  
  return(results)
}

# Function to process all files and save the results
process_all_files_1 <- function(file_list, work_status_list) {
  results <- lapply(names(work_status_list), function(id) {
    file_path <- file_list[grep(id, file_list)]
    if (length(file_path) == 0) {
      message(paste("No file found for ID:", id))
      return(NULL)
    }
    work_status <- work_status_list[[id]]
    message(paste("Processing file:", file_path))
    process_file_1(file_path, work_status)
  })
  
  names(results) <- names(work_status_list)
  return(results)
}





# Splitting the results into work and non-work categories
split_results <- function(results) {
  work_results <- list()
  non_work_results <- list()
  
  for (id in names(results)) {
    if (!is.null(results[[id]]$work)) {
      work_results[[id]] <- results[[id]]$work
    }
    if (!is.null(results[[id]]$non_work)) {
      non_work_results[[id]] <- results[[id]]$non_work
    }
  }
  return(list(work = work_results, non_work = non_work_results))
}





# Function to split data based on work status
split_data_by_work_status_1 <- function(data, work_status) {
  work_days <- which(work_status == 1)
  non_work_days <- which(work_status == 0)
  
  # Handle case where there's only one work day or non-work day
  if (length(work_days) > 0) {
    work_columns <- c("t_values", paste0("MovementAcceleration_g_d", work_days))
  } else {
    work_columns <- c("t_values")
  }
  
  if (length(non_work_days) > 0) {
    non_work_columns <- c("t_values", paste0("MovementAcceleration_g_d", non_work_days))
  } else {
    non_work_columns <- c("t_values")
  }

  
  missing_work_columns <- setdiff(work_columns, colnames(data))
  missing_non_work_columns <- setdiff(non_work_columns, colnames(data))
  
  if (length(missing_work_columns) > 0 && length(work_days) > 0) {
    stop(paste("Some work day columns are missing in the data frame:", paste(missing_work_columns, collapse = ", ")))
  }
  if (length(missing_non_work_columns) > 0 && length(non_work_days) > 0) {
    stop(paste("Some non-work day columns are missing in the data frame:", paste(missing_non_work_columns, collapse = ", ")))
  }
  
  work_data <- data[, work_columns, drop = FALSE]
  non_work_data <- data[, non_work_columns, drop = FALSE]
  
  list(work_data = work_data, non_work_data = non_work_data)
}





# Function to count the number of rows in each data frame
count_rows_in_files <- function(data_list) {
  row_counts <- sapply(data_list, function(data) {
    if (is.null(data)) {
      return(NA)
    } else {
      return(nrow(data))
    }
  })
  return(row_counts)
}





# Function to extract and combine coefs from multiple lists into a single matrix
combine_coefs_into_one <- function(pain_work_lists) {
  
  # Initialize an empty matrix to store the coefficients for all samples
  combined_coefs_matrix <- NULL
  
  # Loop through each list in pain_work_lists
  for (list_index in seq_along(pain_work_lists)) {
    
    # Get the current pain_work list
    current_list <- pain_work_lists[[list_index]]
    
    # Loop through each u-prefixed element in the current list
    for (name in names(current_list)) {
      if (startsWith(name, "u")) {
        # Extract the coefficients for each u-prefixed element
        coefs <- current_list[[name]]$smoothed_mean_fd$fd$coefs
        
        # Combine (bind) these coefficients to the growing matrix
        if (is.null(combined_coefs_matrix)) {
          # Initialize the matrix with the first set of coefficients
          combined_coefs_matrix <- coefs
        } else {
          # Bind the coefficients as additional columns to the existing matrix
          combined_coefs_matrix <- cbind(combined_coefs_matrix, coefs)
        }
      }
    }
  }
  
  return(combined_coefs_matrix)
}






pointwise_sd <- function(beta_fd, XtX_inv_diag, gamma_hat, t_values) {
  # Calculate the pointwise variances
  pointwise_variances <- sapply(t_values, function(t_val) {
    gamma_t_t <- eval.bifd(t_val, t_val, gamma_hat)  # Get the covariance gamma(t, t)
    XtX_inv_diag * gamma_t_t  # Calculate the variance for the time point
  })
  
  # Return the pointwise standard deviations
  sqrt(pointwise_variances)
}





##  compute gamma otimes
compute_gamma_otimes2 <- function(gamma) {
  # Compute gamma^otimes2 using matrix multiplication
  gamma_otimes2 <- gamma %*% gamma
  
  return(gamma_otimes2)
}




# Function to transform difference functions into wide-format data
transform_difference_functions <- function(difference_functions) {
  # Initialize a data frame to collect all difference functions
  all_diff_df <- data.frame(t_values = 1:1440)  # Assuming 1440 time points per day
  
  # Loop through each subject ID to format the difference functions
  for (id in names(difference_functions)) {
    # Get the difference function for the current subject
    diff_function <- difference_functions[[id]]
    
    # Create a data frame for the current subject's difference function
    subject_diff_df <- data.frame(
      t_values = 1:1440,  # Time values (assuming one-minute intervals)
      Difference = diff_function
    )
    
    # Pivot the data to create wide format with one column per subject
    subject_diff_wide <- subject_diff_df %>%
      pivot_wider(names_from = t_values, values_from = Difference, names_prefix = paste0(id, "_"))
    
    # Merge the current subject's data with the main data frame
    all_diff_df <- cbind(all_diff_df, subject_diff_wide)
  }
  
  # Return the final wide-format data frame
  return(all_diff_df)
}




### testing
perform_tests <- function(beta_list, residual_values, XtX_inv_diag, n_basis, n, p_interaction) {
  
  # Calculate covariance matrix for residuals
  V_hat <- cov(residual_values)
  
  # Initialize result lists
  K_values <- numeric(length(beta_list))
  p_values_K <- numeric(length(beta_list))
  significant_K <- logical(length(beta_list))
  
  p_values_WS <- numeric(length(beta_list))
  significant_L2 <- logical(length(beta_list))
  
  p_values_F <- numeric(length(beta_list))
  significant_F <- logical(length(beta_list))
  
  # Perform tests for each beta function
  for (i in seq_along(beta_list)) {
    beta_hat <- beta_list[[i]]
    
    # K Test
    K_value <- c(c(beta_hat) %*% solve(V_hat) / XtX_inv_diag[i + 1]) %*% c(beta_hat)
    p_value_K <- 1 - pchisq(K_value, df = n_basis)
    significant_K[i] <- K_value > qchisq(0.95, df = n_basis)
    K_values[i] <- K_value
    p_values_K[i] <- p_value_K
    
    # L2 Test (Welch-Satterthwaite)
    Cov_beta <- V_hat * XtX_inv_diag[i + 1]
    gamma_otimes2_beta <- compute_gamma_otimes2(Cov_beta)
    
    trace_gamma_beta <- sum(diag(Cov_beta))
    trace_gamma_otimes2_beta <- sum(diag(gamma_otimes2_beta))
    
    beta_hat_k <- trace_gamma_otimes2_beta / trace_gamma_beta
    k_hat_k <- (trace_gamma_beta^2) / trace_gamma_otimes2_beta
    
    Tn_K_actual <- t(beta_hat) %*% beta_hat
    critical_value_2 <- beta_hat_k * qchisq(0.95, df = k_hat_k)
    
    p_value_L2 <- 1 - pchisq(Tn_K_actual / beta_hat_k, df = k_hat_k)
    significant_L2[i] <- Tn_K_actual > critical_value_2
    
    p_values_WS[i] <- p_value_L2
    
    # F Test
    Fn_value <- Tn_K_actual / trace_gamma_beta
    
    df1_k <- k_hat_k
    df2_k <- (n - p_interaction - 1) * k_hat_k
    
    critical_value_3 <- qf(0.95, df1 = df1_k, df2 = df2_k)
    significant_F[i] <- Fn_value > critical_value_3
    
    p_value_F <- 1 - pf(Fn_value, df1 = df1_k, df2 = df2_k)
    p_values_F[i] <- p_value_F
  }
  
  # Return results
  return(list(
    K_test_p_values = p_values_K,
    K_test_significant = significant_K,
    L2_test_p_values = p_values_WS,
    L2_test_significant = significant_L2,
    F_test_p_values = p_values_F,
    F_test_significant = significant_F
  ))
}





create_covariance_structure <- function(t_values, V) {
  k <- nrow(V)  # Number of basis functions
  basis <- create_basis_functions(k, "fourier", rangeval = range(t_values))  # Create the basis functions
  Phi <- eval.basis(t_values, basis)  # Evaluate the basis functions at the time points
  
  # Initialize the covariance matrix
  Sigma <- matrix(0, nrow = length(t_values), ncol = length(t_values))
  
  # Calculate the covariance function as per the formula
  for (h in 1:k) {
    for (h_prime in 1:k) {
      Sigma <- Sigma + V[h, h_prime] * outer(Phi[, h], Phi[, h_prime])
    }
  }
  
  # Force Sigma to be positive semi-definite
  Sigma <- force_positive_semidefinite(Sigma)
  
  return(Sigma)
}





### create basis functions
create_basis_functions <- function(k = 6, basis_type = "fourier", rangeval = c(0, 5)) {
  if (basis_type == "fourier") {
    basis <- create.fourier.basis(rangeval = rangeval, nbasis = k)
  }
  return(basis)
}




force_positive_semidefinite <- function(matrix) {
  eigen_decomp <- eigen(matrix)
  eigen_decomp$values <- pmax(eigen_decomp$values, 0.00001)  # Set any negative eigenvalues to zero
  return(eigen_decomp$vectors %*% diag(eigen_decomp$values) %*% t(eigen_decomp$vectors))
}
##without force you get eigenvalues like -7.888899e-18









calculate_steps_per_day <- function(file_path, sep = ";", dec = ",", xlsx = TRUE) {
  # Read and transform data
  data_wide <- read_and_transform_data(file_path, variable = "StepCount_steps")
  
  # Summing steps per day
  step_sums <- colSums(data_wide[,-1], na.rm = TRUE) # Exclude the time column (first column) from summation
  
  return(step_sums)
}

# Function to process all files and calculate steps for each person
process_all_files_steps <- function(file_list, work_status_list) {
  results <- lapply(names(work_status_list), function(id) {
    file_path <- file_list[grep(id, file_list)]
    if (length(file_path) == 0) {
      message(paste("No file found for ID:", id))
      return(NULL)
    }
    
    message(paste("Processing file:", file_path))
    total_steps <- calculate_steps_per_day(file_path)
    
    # Combine results with person ID
    return(data.frame(ID = id, t(total_steps)))
  })
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results)
  
  return(results_df)
}





calculate_mean_work_nonwork_steps <- function(results_steps, work_status_list) {
  results <- lapply(1:nrow(results_steps), function(i) {
    # Extract the ID for the current person
    id <- results_steps$ID[i]
    
    # Extract the work status for the given ID
    work_status <- work_status_list[[id]]
    
    # Create logical indices for workdays and non-work days based on work_status
    work_indices <- which(work_status == 1)
    non_work_indices <- which(work_status == 0)
    
    # Extract the columns for workdays and non-work days (skip the ID column)
    work_columns <- results_steps[i, work_indices + 1]  # Adjust for ID column
    non_work_columns <- results_steps[i, non_work_indices + 1]
    
    # Initialize mean values as NA
    work_mean <- NA
    non_work_mean <- NA
    
    # Check if there are any work day columns and handle vector vs matrix case
    if (length(work_indices) > 0 && !all(is.na(work_columns))) {
      if (is.vector(work_columns) || ncol(work_columns) == 1) {
        work_mean <- mean(work_columns, na.rm = TRUE)
      } else {
        work_mean <- rowMeans(work_columns, na.rm = TRUE)
      }
    }
    
    # Check if there are any non-work day columns and handle vector vs matrix case
    if (length(non_work_indices) > 0 && !all(is.na(non_work_columns))) {
      if (is.vector(non_work_columns) || ncol(non_work_columns) == 1) {
        non_work_mean <- mean(non_work_columns, na.rm = TRUE)
      } else {
        non_work_mean <- rowMeans(non_work_columns, na.rm = TRUE)
      }
    }
    
    return(data.frame(ID = id, work_mean = work_mean, non_work_mean = non_work_mean))
  })
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results)
  
  return(results_df)
}










################################################################################
#############--------PLOTS----------############################################
################################################################################


# Define a custom theme with no background grid
custom_theme <- theme(
  panel.background = element_blank(),  # Remove background
  panel.grid.major = element_line(size = 0.5),
  panel.grid.minor = element_line(size = 0.25),
  axis.line = element_line(colour = "black", size = 0.5),  # Axis lines
  axis.text.x = element_text(size = 22),  # Adjusted size for x-axis text
  axis.text.y = element_text(size = 22),  # Adjusted size for y-axis text
  axis.title.x = element_text(size = 22),  # Adjusted size for x-axis title
  axis.title.y = element_text(size = 22),  # Adjusted size for y-axis title
  plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),  # Adjusted size, bold and centered plot title
  legend.title = element_blank(),  # No legend title
  legend.text = element_text(size = 20),  # Adjusted size for legend text
  legend.background = element_blank(),  # No background for the legend
  legend.key = element_blank(),  # No boxes around the legend keys
  legend.position = "bottom",  # Position legend to the right
  strip.text = element_text(face = "bold", size = 22),  # Facet strip text size
  strip.background = element_blank(),  # Remove background box
  panel.spacing = unit(0.5, "lines"),  # Space between panels
  panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Single line border around each panel
)



# Combined Plot Function without shaded areas
plot_combined_functions_with_mean <- function(smoothed_data_work, smoothed_data_non_work) {
  # Evaluate smoothed data
  eval_points <- seq(0, 1440, length.out = 1440)
  smoothed_values_work <- eval.fd(eval_points, smoothed_data_work$fd)
  smoothed_values_non_work <- eval.fd(eval_points, smoothed_data_non_work$fd)
  
  # Convert time from minutes to hours
  time_hours <- eval_points / 60
  
  # Create data frames for ggplot
  work_df <- as.data.frame(smoothed_values_work) %>%
    mutate(Time = time_hours, Group = "Work Days")
  
  non_work_df <- as.data.frame(smoothed_values_non_work) %>%
    mutate(Time = time_hours, Group = "Non-Work Days")
  
  combined_df <- bind_rows(work_df, non_work_df)
  
  # Pivot data for ggplot
  combined_long_df <- combined_df %>%
    pivot_longer(cols = -c(Time, Group), names_to = "Function", values_to = "Value")
  
  # Plot
  ggplot(combined_long_df, aes(x = Time, y = Value, color = Group)) +
    geom_line(aes(group = interaction(Group, Function)), alpha = 0.3, size = 0.5) +  # Thin lines for each day
    stat_summary(fun = mean, geom = "line", aes(group = Group), size = 1.5, linetype = "solid") +  # Thicker line for mean
    labs(
      title = "Combined Smoothed Movement Acceleration",
      y = "Movement Acceleration (g)",
      x = "Time (hours)"
    ) +
    scale_x_continuous(breaks = seq(0, 24, by = 2)) +
    scale_color_manual(values = c("Work Days" = "blue", "Non-Work Days" = "red")) +
    custom_theme +
    guides(color = guide_legend(override.aes = list(size = 1.5, linetype = "solid")))
}



plot_boxplots_comparison <- function(group1_data, group2_data, group1_label = "Group 1", group2_label = "Group 2", title = "Comparison", y_label = "Steps") {
  
  # Combine the data into a data frame for plotting, and divide by 1000 for 'thousands' steps
  data <- data.frame(
    Steps = c(group1_data, group2_data) / 1000,  # Steps in thousands
    Group = factor(c(rep(group1_label, length(group1_data)), rep(group2_label, length(group2_data))))
  )
  
  # Filter out NA values
  data <- data[!is.na(data$Steps), ]
  
  # Calculate statistics for group 1
  group1_count <- sum(!is.na(group1_data))
  group1_mean <- mean(group1_data, na.rm = TRUE) / 1000
  group1_median <- median(group1_data, na.rm = TRUE) / 1000
  
  # Calculate statistics for group 2
  group2_count <- sum(!is.na(group2_data))
  group2_mean <- mean(group2_data, na.rm = TRUE) / 1000
  group2_median <- median(group2_data, na.rm = TRUE) / 1000
  
  # Create the boxplot
  boxplot <- ggplot(data, aes(x = Group, y = Steps, fill = Group)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
    labs(title = title, y = y_label, x = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 16),  # Larger x-axis text
      axis.text.y = element_text(size = 16),  # Larger y-axis text
      axis.title.x = element_text(size = 16),  # Larger x-axis title
      axis.title.y = element_text(size = 16),  # Larger y-axis title
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),  # Larger and centered title
      legend.position = "none"  # Remove the legend
    ) +
    scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "k"))  # Y-axis in thousands
  
  # Return the statistics and plot in a list
  return(list(
    group1_stats = list(
      label = group1_label,
      count = group1_count,
      mean = group1_mean,
      median = group1_median
    ),
    group2_stats = list(
      label = group2_label,
      count = group2_count,
      mean = group2_mean,
      median = group2_median
    ),
    plot = boxplot
  ))
}



# Plotting function for covariance matrices
plot_covariance_matrix <- function(Sigma, title, x = "t", y = "t'") {
  Sigma_melted <- melt(Sigma)
  
  plot <- ggplot(Sigma_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                         limit = c(min(Sigma_melted$value), max(Sigma_melted$value)), 
                         breaks = c(min(Sigma_melted$value), max(Sigma_melted$value)),  # Display only min and max in legend
                         labels = round(c(min(Sigma_melted$value), max(Sigma_melted$value)), 2)) +  # Format labels
    scale_x_continuous(limits = c(0, 1440), expand = c(0, 0)) +  # Set x-axis from 0 to 1440
    scale_y_continuous(limits = c(0, 1440), expand = c(0, 0)) +  # Set y-axis from 0 to 1440
    labs(title = title, x = x, y = y) +
    custom_theme +
    theme(
      panel.border = element_blank(),  # Remove the panel border
      axis.line = element_blank()  # Remove axis lines to avoid the box
    )
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}



# Function to determine the optimal number of basis functions using the Scree Plot Method
determine_optimal_nbasis <- function(sim_data, basis_type = "fourier", max_basis = 100) {
  # Extract time values
  t_values <- sim_data$t_values
  
  # Initialize list to store RSS values and plots for each subject
  rss_values_list <- list()
  plot_list <- list()
  optimal_basis_numbers <- integer(ncol(sim_data) - 1)
  
  # Loop over each observation/subject (each column except the first one which is t_values)
  for (i in 2:ncol(sim_data)) {
    single_subject_data <- sim_data[[i]]  # Extract data for one observation
    
    # Initialize vector to store RSS values
    rss_values <- numeric(max_basis - 4 + 1)
    
    # Loop over a range of basis numbers to calculate RSS
    for (nbasis in 4:max_basis) {
      # Create basis according to the specified type
      if (basis_type == "fourier") {
        basis <- create.fourier.basis(c(min(t_values), max(t_values)), nbasis)
      } else if (basis_type == "bspline") {
        basis <- create.bspline.basis(c(min(t_values), max(t_values)), nbasis)
      } else if (basis_type == "constant") {
        basis <- create.constant.basis(c(min(t_values), max(t_values)), nbasis)
      } else if (basis_type == "monomial") {
        basis <- create.monomial.basis(c(min(t_values), max(t_values)), nbasis)
      } else {
        stop("Invalid basis_type. Please choose 'fourier', 'bspline', or 'monomial'.")
      }
      
      # Fit the model and calculate residual sum of squares
      smooth_result <- smooth.basis(t_values, single_subject_data, basis)
      residuals <- single_subject_data - eval.fd(t_values, smooth_result$fd)
      rss <- sum(residuals^2)
      rss_values[nbasis - 4 + 1] <- rss
    }
    
    # Store the RSS values for the current subject
    rss_values_list[[i - 1]] <- rss_values
    
    # Create and store the Scree Plot for each subject
    nbasis_values <- seq(4, max_basis)
    rss_data <- data.frame(nbasis = nbasis_values, rss = rss_values)
    rss_plot <- ggplot(rss_data, aes(x = nbasis, y = rss)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Scree Plot"),  # Updated to include subject number
           x = "Number of Basis Functions",
           y = "Residual Sum of Squares (RSS)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))  # Centering the title
    
    plot_list[[i - 1]] <- rss_plot  # Store the plot in the list
  }
  
  return(plot_list)  # Return the list of plots
}




export_plot_with_ci <- function(filename, beta_fd, se_beta, title, color, eval_grid) {
  png(filename, width = 800, height = 800, res = 130)
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot_with_ci(beta_fd, se_beta, title, color, eval_grid)
  dev.off()
}



export_plot_with_ci_age <- function(filename, beta_fd, se_beta, title, color, eval_grid) {
  png(filename, width = 800, height = 800, res = 130)
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot_with_ci_age(beta_fd, se_beta, title, color, eval_grid)
  dev.off()
}



# Plotting function for V matrices with adjusted axis labels
plot_V_matrix <- function(V_matrix, title, x = expression(epsilon[h]), y = expression(epsilon[h])) {
  V_melted <- melt(V_matrix)
  
  # Check if the variable is numeric or discrete
  is_numeric <- is.numeric(V_melted$Var1)
  
  # Adjust scale accordingly
  plot <- ggplot(V_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                         limit = c(min(V_melted$value), max(V_melted$value))) +
    labs(title = title, x = x, y = y) +
    labs(title = title, x = x, y = y) +
    ### general theme
    # Define a custom theme with no background grid
    theme(
      panel.background = element_blank(),  # Remove background
      panel.grid.major = element_line(size = 0.5),
      panel.grid.minor = element_line(size = 0.25),
      axis.line = element_blank(),  # Remove axis lines to avoid the box
      axis.text.x = element_text(size = 15),  # Adjusted size for x-axis text
      axis.text.y = element_text(size = 15),  # Adjusted size for y-axis text
      axis.title.x = element_text(size = 22),  # Adjusted size for x-axis title
      axis.title.y = element_text(size = 22),  # Adjusted size for y-axis title
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),  # Adjusted size, bold and centered plot title
      legend.title = element_blank(),  # No legend title
      legend.text = element_text(size = 20),  # Adjusted size for legend text
      legend.background = element_blank(),  # No background for the legend
      legend.key = element_blank(),  # No boxes around the legend keys
      legend.position = "right",  # Position legend to the right
      strip.text = element_text(face = "bold", size = 22),  # Facet strip text size
      strip.background = element_blank(),  # Remove background box
      panel.spacing = unit(0.5, "lines"),  # Space between panels
      panel.border = element_blank()  # Remove the panel border
    )
  if (!is_numeric) {
    plot <- plot +
      scale_x_discrete(breaks = seq(0, 34, by = 2), labels = as.character(seq(0, 34, by = 2))) +  # Change X-axis labels
      scale_y_discrete(breaks = seq(0, 34, by = 2), labels = as.character(seq(0, 34, by = 2)))    # Change Y-axis labels
  }
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}









# Plotting function for covariance matrices
plot_covariance_matrix <- function(Sigma, title, x = "" , y = "") {
  Sigma_melted <- melt(Sigma)
  plot <- ggplot(Sigma_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                         limit = c(min(Sigma_melted$value), max(Sigma_melted$value))) +
    labs(title = title, x = x, y = y) +
    ### general theme
    # Define a custom theme with no background grid
    theme(
      panel.background = element_blank(),  # Remove background
      panel.grid.major = element_line(size = 0.5),
      panel.grid.minor = element_line(size = 0.25),
      axis.line = element_blank(),  # Remove axis lines to avoid the box
      axis.text.x = element_blank(),  # Adjusted size for x-axis text
      axis.text.y = element_blank(),  # Adjusted size for y-axis text
      axis.title.x = element_text(size = 22),  # Adjusted size for x-axis title
      axis.title.y = element_text(size = 22),  # Adjusted size for y-axis title
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),  # Adjusted size, bold and centered plot title
      legend.title = element_blank(),  # No legend title
      legend.text = element_text(size = 20),  # Adjusted size for legend text
      legend.background = element_blank(),  # No background for the legend
      legend.key = element_blank(),  # No boxes around the legend keys
      legend.position = "right",  # Position legend to the right
      strip.text = element_text(face = "bold", size = 22),  # Facet strip text size
      strip.background = element_blank(),  # Remove background box
      panel.spacing = unit(0.5, "lines"),  # Space between panels
      panel.border = element_blank()  # Remove the panel border
    )
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}





pointwise_sd_1 <- function(beta_fd, XtX_inv_diag, gamma_hat, t_values) {
  # Get the number of time points
  n_time_points <- length(t_values)
  
  # Create a matrix to store pointwise variances
  pointwise_variances <- numeric(n_time_points)
  
  for (i in 1:n_time_points) {
    t_val <- t_values[i]
    # Evaluate the covariance function at all pairs (t_val, t_val)
    gamma_t_t <- eval.bifd(t_val, t_val, gamma_hat)
    
    # Consider interactions with neighboring points (smoothing across time points)
    for (j in 1:n_time_points) {
      if (i != j) {
        t_val_other <- t_values[j]
        gamma_t_other <- eval.bifd(t_val, t_val_other, gamma_hat)
        gamma_t_t <- gamma_t_t + 2 * gamma_t_other
      }
    }
    
    # Calculate the variance for the time point, scaling with XtX_inv_diag
    pointwise_variances[i] <- XtX_inv_diag * gamma_t_t
  }
  
  # Return the pointwise standard deviations
  sqrt(pointwise_variances)
}



plot_with_ci <- function(beta_fd, se_beta, title, color = "blue", eval_grid) {
  # Set up the plotting area with appropriate margins
  par(mar = c(4, 4, 4, 2) + 1, mgp = c(1.5, 0.5, 0))  # Adjust mgp to bring labels closer
  
  # Determine the beta values and their confidence intervals
  beta_hat <- as.numeric(eval.fd(beta_fd, eval_grid))
  ci_lower <- beta_hat - 1.96 * se_beta
  ci_upper <- beta_hat + 1.96 * se_beta
  
  # Define y-axis range for consistency
  y_range <- range(c(ci_lower, ci_upper), na.rm = TRUE)
  
  # Plot the beta function with confidence intervals
  plot(eval_grid, beta_hat, type = "l", col = color, lwd = 2, ylim = y_range, main = title,
       xlab = "Time (hours)", ylab = "MET", xlim = c(0, 1440), xaxt = "n", cex.axis = 1, cex.lab =1 , cex.main = 1)
  
  # Add confidence band
  polygon(c(eval_grid, rev(eval_grid)), c(ci_lower, rev(ci_upper)), col = adjustcolor(color, alpha.f = 0.2), border = NA)
  lines(eval_grid, beta_hat, col = color, lwd = 3)
  
  # Add a horizontal reference line at zero
  abline(h = 0, col = "red", lty = 2)
  
  # Add custom x-axis labels for every hour
  axis(1, at = seq(0, 1440, by = 120), labels = seq(0, 24, by = 2), cex.axis = 1, las = 1)
  
  # Add the legend below the plot, adjusted for proper display
  # Add the legend below the plot
  #legend("bottom", legend = c("Estimated Beta Function", "Confidence Band", "Zero Line"),
  #       col = c(color, "gray70", "red"),
  #       lwd = c(3, 8, 1), lty = c(1, 1, 2), bty = "o", inset = c(0, -0.2), xpd = TRUE, cex = 1.4, horiz = TRUE)
}




plot_with_ci_age <- function(beta_fd, se_beta, title, color = "blue", eval_grid) {
  # Set up the plotting area with appropriate margins
  par(mar = c(4, 4, 4, 2) + 1, mgp = c(1.5, 0.5, 0))  # Adjust mgp to bring labels closer
  
  # Determine the beta values and their confidence intervals
  beta_hat <- as.numeric(eval.fd(beta_fd, eval_grid))
  ci_lower <- beta_hat - 1.96 * se_beta
  ci_upper <- beta_hat + 1.96 * se_beta
  
  # Define y-axis range for consistency
  y_range <- range(c(ci_lower, ci_upper), na.rm = TRUE)
  
  # Plot the beta function with confidence intervals
  plot(eval_grid, beta_hat, type = "l", col = color, lwd = 2, ylim = c(-0.01,0.025), main = title,
       xlab = "Time (hours)", ylab = "MET", xlim = c(0, 1440), xaxt = "n", cex.axis = 1, cex.lab =1 , cex.main = 1)
  
  # Add confidence band
  polygon(c(eval_grid, rev(eval_grid)), c(ci_lower, rev(ci_upper)), col = adjustcolor(color, alpha.f = 0.2), border = NA)
  lines(eval_grid, beta_hat, col = color, lwd = 3)
  
  # Add a horizontal reference line at zero
  abline(h = 0, col = "red", lty = 2)
  
  # Add custom x-axis labels for every hour
  axis(1, at = seq(0, 1440, by = 120), labels = seq(0, 24, by = 2), cex.axis = 1, las = 1)
  
  # Add the legend below the plot, adjusted for proper display
  # Add the legend below the plot
  #legend("bottom", legend = c("Estimated Beta Function", "Confidence Band", "Zero Line"),
  #       col = c(color, "gray70", "red"),
  #       lwd = c(3, 8, 1), lty = c(1, 1, 2), bty = "o", inset = c(0, -0.2), xpd = TRUE, cex = 1.4, horiz = TRUE)
}



# Define the function to plot Pain vs. No Pain for a given group, with customizable legend labels
plot_pain_vs_no_pain <- function(fd_pain, fd_nopain, mean_fd_pain, mean_fd_nopain, title, 
                                 legend_labels = c("a", "a", "a", "a"), y_range = range(0.98, 2.7)) {
  # Set up the plotting area with appropriate margins and adjusted mgp
  par(mar = c(6, 4, 4, 2) + 1, mgp = c(1.5, 0.5, 0))  # Adjust mgp to bring labels closer
  
  # Determine the shared y-axis range for consistency
  y_range <- y_range
  
  # Plot Pain vs. No Pain
  plot(fd_pain, main = title, col = adjustcolor("blue", alpha.f = 0.1), xlab = "Time (hours)", ylab = "MET",
       xlim = c(0, 1440), ylim = y_range, xaxt = "n", cex.axis = 1)  # Suppress default x-axis and reduce label size
  lines(mean_fd_pain, col = "blue", lwd = 3)
  lines(fd_nopain, col = adjustcolor("red", alpha.f = 0.1))
  lines(mean_fd_nopain, col = "red", lwd = 3)
  
  # Add custom x-axis labels for every hour
  axis(1, at = seq(0, 1440, by = 120), labels = seq(0, 24, by = 2), cex.axis = 1, las = 1)  # las = 1 for horizontal labels
  
  # Add the legend below the plot inside the plot area, in one row, with customizable legend labels
  legend("bottom", legend = legend_labels,
         col = c("blue", "red", adjustcolor("blue", alpha.f = 0.2), adjustcolor("red", alpha.f = 0.2)),
         lwd = c(3, 3, 1, 1), ncol = 4, bty = "n", inset = c(0, -1.2), xpd = TRUE, cex = 1)
}



# Function to plot smoothed mean functions for different comparisons
plot_smoothed_mean_comparisons <- function(overall_mean_results_list, titles, colors, plot_title) {
  # Initialize an empty list to collect data for plotting
  combined_data <- data.frame()
  
  # Loop over each result and extract smoothed mean functions
  for (i in seq_along(overall_mean_results_list)) {
    overall_mean_result <- overall_mean_results_list[[i]]
    title <- titles[i]
    
    time_hours <- 1:1440 / 60  # Assuming time is in minutes, convert to hours
    
    # Extract smoothed mean function data
    smoothed_data <- data.frame(
      time = time_hours,
      value = eval.fd(1:1440, overall_mean_result$smoothed_overall_mean_fd$fd),
      group = title
    )
    
    # Combine smoothed data
    combined_data <- rbind(combined_data, smoothed_data)
  }
  
  # Plot using ggplot2
  ggplot(combined_data, aes(x = time, y = value, color = group)) +
    geom_line(size = 1) +
    labs(title = plot_title,
         x = "Time (hours)",
         y = "Movement Acceleration (g)") +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = seq(0, 24, by = 2)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 19, face = "bold"),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 17),
      legend.position = "bottom"  # Position legend underneath the plot
    )
}






### . Plot Raw and Smoothed Data
plot_raw_and_smoothed_data <- function(raw_data, smoothed_data, day, title = "Raw and smoothed data") {
  smoothed_day <- eval.fd(seq(0, 1440, length.out = 1440), smoothed_data$fd)[, day]
  
  plot_data <- data.frame(
    time_point = seq(0, 1440, length.out = 1440) / 60,  # Convert minutes to hours
    raw_data = raw_data,
    smoothed_data = smoothed_day
  )
  
  ggplot(plot_data, aes(x = time_point)) +
    geom_line(aes(y = raw_data, color = "Raw Data"), size = 0.3) +
    geom_line(aes(y = smoothed_data, color = "Smoothed Data"), linetype = "solid", size = 2) +
    labs(title = title,
         y = "Movement Acceleration (g)",
         x = "Time (hours)") +
    scale_x_continuous(breaks = seq(0, 24, by = 2)) +  # Adjust x-axis breaks to show hours
    scale_color_manual(values = c("Raw Data" = "blue", "Smoothed Data" = "red")) +
    custom_theme +  # Apply the custom theme here
    guides(color = guide_legend(override.aes = list(size = 2, linetype = "solid")))
}


### Plotting Functions for Log-Transformed Data
plot_log_data <- function(raw_data_log, smoothed_data_log, day, title = "Log-Transformed Movement Acceleration Data") {
  smoothed_day_log <- eval.fd(seq(0, 1440, length.out = 1440), smoothed_data_log$fd)[, day]
  
  plot_data_log <- data.frame(
    time_point = seq(0, 1440, length.out = 1440),
    log_raw_data = raw_data_log,
    log_smoothed_data = smoothed_day_log
  )
  
  ggplot(plot_data_log, aes(x = time_point)) +
    geom_line(aes(y = log_raw_data, color = "Log Raw Data"), size = 0.3) +
    geom_line(aes(y = log_smoothed_data, color = "Log Smoothed Data"), linetype = "solid", size = 1) +
    labs(title = title,
         y = "Log Movement Acceleration (g)",
         x = "Time (hours)") +
    scale_x_continuous(breaks = seq(0, 24, by = 2), labels = seq(0, 24, by = 2)) +  # Adjust x-axis breaks to show hours
    scale_color_manual(values = c("Log Raw Data" = "blue", "Log Smoothed Data" = "red")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 19, face = "bold"),  # Increase title size
      axis.title = element_text(size = 17),  # Increase axis title size
      axis.text = element_text(size = 15),  # Increase axis text size
      legend.text = element_text(size = 15),  # Increase legend text size
      legend.title = element_text(size = 17)  # Increase legend title size
    ) +
    guides(color = guide_legend(override.aes = list(size = 1, linetype = "solid")))
}



###  Create Mean Plot
calculate_mean_function <- function(data) {
  mean_function <- rowMeans(data[, -1], na.rm = TRUE)
  return(mean_function)
}

plot_all_functions_with_mean <- function(smoothed_data, title = "All Functions with Mean Function") {
  mean_fd <- mean.fd(smoothed_data$fd)
  
  eval_points <- seq(0, 1440, length.out = 1440)
  smoothed_values <- eval.fd(eval_points, smoothed_data$fd)
  ylim_max <- 1.1 * max(smoothed_values, na.rm = TRUE)
  
  # Convert time from minutes to hours
  time_hours <- eval_points / 60
  
  # Plot the smoothed data functions with custom x-axis
  plot(time_hours, smoothed_values[,1], type = "n", col = "grey", 
       ylab = "Movement Acceleration (g)", xlab = "", ylim = c(0, ylim_max), xaxt = 'n', 
       main = title, cex.main = 1.7, cex.lab = 1.5, cex.axis = 1.3)
  matlines(time_hours, smoothed_values, col = rgb(0.8, 0.8, 0.8, 0.5), lty = 1)
  
  # Add the mean function line
  lines(time_hours, eval.fd(eval_points, mean_fd), col = "red", lwd = 2)
  
  # Customize x-axis with timepoints every two hours
  axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2), cex.axis = 1.5)
  mtext("Time (hours)", side = 1, line = 2.5, cex = 1.7)
  
  # Add a legend with customized text size
  legend("topright", legend = c("Mean Function"), col = c("red"), lwd = 2, cex = 1.5)
}




plot_all_functions_with_mean_and_smoothed_mean <- function(data, mean_function, smoothed_mean_fd, title = "Raw, Mean, and Smoothed Mean Functions") {
  num_raw_series <- ncol(data) - 1
  
  plot_data <- data.frame(
    time_point = rep(data$t_values, num_raw_series) / 60,  # Convert to hours
    value = as.vector(as.matrix(data[, -1])),
    group = rep(names(data)[-1], each = nrow(data))
  )
  
  mean_plot_data <- data.frame(
    time_point = data$t_values / 60,  # Convert to hours
    value = mean_function,
    group = "Mean Function"
  )
  
  smoothed_mean_values <- eval.fd(data$t_values, smoothed_mean_fd$fd)
  smoothed_mean_plot_data <- data.frame(
    time_point = data$t_values / 60,  # Convert to hours
    value = smoothed_mean_values,
    group = "Smoothed Mean Function"
  )
  
  combined_plot_data <- rbind(plot_data, mean_plot_data, smoothed_mean_plot_data)
  
  ylim_max <- max(mean_function, na.rm = TRUE) * 1.1
  
  ggplot(combined_plot_data, aes(x = time_point, y = value, group = group, color = group)) +
    geom_line(data = subset(combined_plot_data, group %in% names(data)[-1]), alpha = 0.5, aes(color = "Raw Data")) +
    geom_line(data = mean_plot_data, aes(x = time_point, y = value, color = "Mean Function"), size = 1.3) +
    geom_line(data = smoothed_mean_plot_data, aes(x = time_point, y = value, color = "Smoothed Mean Function"), size = 1.3) +
    labs(title = title, y = "Movement Acceleration (g)", x = "Time (hours)") +
    scale_x_continuous(breaks = seq(0, 24, by = 2)) +  # Adjust x-axis breaks to show hours
    scale_color_manual(values = c("Raw Data" = "gray", "Mean Function" = "blue", "Smoothed Mean Function" = "red"),
                       breaks = c("Raw Data", "Mean Function", "Smoothed Mean Function"),
                       labels = c("Raw Data", "Mean Function", "Smoothed Mean Function")) +
    guides(color = guide_legend(override.aes = list(size = c(0.7, 1.3, 1.3), alpha = c(1, 1, 1)))) +  # Adjust line size in legend
    ylim(0, ylim_max) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 19, face = "bold"),  # Increase title size
      axis.title = element_text(size = 17),  # Increase axis title size
      axis.text = element_text(size = 15),  # Increase axis text size
      legend.text = element_text(size = 15),  # Increase legend text size
      legend.title = element_text(size = 17)  # Increase legend title size
    )
}





plot_mean_and_smoothed_mean <- function(result, title) {
  # Convert time from minutes to hours
  time_hours <- 1:1440 / 60
  
  # Plot the mean function
  plot(time_hours, result$mean_function_raw, type = "l", col = "blue", lwd = 2, 
       ylab = "Movement Acceleration (g)", xlab = "", 
       main = title, xaxt = 'n', cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.3)
  
  # Customize x-axis with timepoints every two hours
  axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2), cex.axis = 1.3)
  mtext("Time (hours)", side = 1, line = 2.5, cex = 1.3)
  
  # Add the smoothed mean function line
  smoothed_values <- eval.fd(1:1440, result$smoothed_mean_fd$fd)
  lines(time_hours, smoothed_values, col = "red", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Mean Function", "Smoothed Mean Function"), col = c("blue", "red"), lwd = 2, cex = 1.2)
}




plot_all_raw_and_overall_mean_functions_old <- function(results, overall_mean_result, title = "Raw Mean Functions and Overall Mean Function") {
  # Combine all raw mean functions into a data frame for plotting
  all_raw_means_df <- do.call(cbind, lapply(results, function(res) res$mean_function_raw))
  
  # Convert time from minutes to hours
  time_hours <- 1:1440 / 60
  
  # Plot all raw mean functions
  matplot(time_hours, all_raw_means_df, type = "l", col = rgb(0.8, 0.8, 0.8, 0.5), lty = 1, 
          ylab = "Movement Acceleration (g)", xlab = "", 
          main = title, 
          ylim = c(0, 1), xaxt = 'n', cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.3)
  
  # Customize x-axis with timepoints every two hours
  axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2), cex.axis = 1.3)
  mtext("Time (hours)", side = 1, line = 2.5, cex = 1.7)
  
  # Add the overall raw mean function
  lines(time_hours, overall_mean_result$overall_mean_function, col = "blue", lwd = 2)
  
  # Add the smoothed overall mean function
  lines(time_hours, eval.fd(1:1440, overall_mean_result$smoothed_overall_mean_fd$fd), col = "red", lwd = 2)
  
  # Add a legend with customized text size
  legend("topright", legend = c("Individual Mean Functions", "Overall Mean Function", "Smoothed Overall Mean Function"), 
         col = c(rgb(0.8, 0.8, 0.8, 0.5), "blue", "red"), lwd = 2, cex = 1.3)}




plot_all_raw_and_overall_mean_functions <- function(results, title = "Raw Mean Functions and Overall Mean Function") {
  # Combine all raw mean functions into a data frame for plotting
  all_raw_means_df <- do.call(cbind, lapply(results, function(res) res$mean_function_raw))
  
  # Convert time from minutes to hours
  time_hours <- 1:1440 / 60
  
  # Set up the plot with increased margins for the legend and axes labels
  par(mar = c(7, 7, 4, 2) + 0.1)  # Increase bottom and left margins for spacing
  plot(time_hours, all_raw_means_df[, 1], type = "n", col = "grey", lty = 1, 
       ylab = "Movement Acceleration (g)", xlab = "Time (hours)", 
       main = title, 
       ylim = c(0, 4), xaxt = 'n', cex.main = 2, cex.lab = 2, cex.axis = 2,
       mgp = c(4, 1, 0))  # Move axis labels away from the axis
  
  # Plot all raw mean functions as thin lines
  matlines(time_hours, all_raw_means_df, col = "grey", lty = 1, lwd = 0.5)  # Thinner lines for individual functions
  
  # Customize x-axis with timepoints every two hours
  axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2), cex.axis = 1.6, mgp = c(3, 1, 0))  # Move x-axis labels further from axis
  
  
  # Add a legend with customized text size underneath the "Time (hours)" label
  legend("topright", legend = c("Individual raw mean functions"), 
         col = c("grey"), lwd = c(0.5), cex = 1.6, horiz = FALSE, 
         bty = "o", inset = c(0,0.01), x.intersp = 0.5, y.intersp = 0.5, seg.len = 2)
  
  # Reset margins to default
  par(mar = c(5, 4, 4, 2) + 0.1)
}





# Function to plot smoothed mean functions for pain and no pain with adjustable linewidth
plot_combined_smoothed_mean_functions <- function(pain_smoothed_mean_fd, no_pain_smoothed_mean_fd, line_width = 1.5) {
  # Combine the smoothed mean functions into one data frame
  combined_means <- rbind(
    combine_smoothed_mean_functions(pain_smoothed_mean_fd, "Pain"),
    combine_smoothed_mean_functions(no_pain_smoothed_mean_fd, "No Pain")
  )
  
  combined_means$t_values <- combined_means$t_values / 60
  
  # Plot the combined smoothed mean functions
  ggplot(combined_means, aes(x = t_values, y = mean_function, color = label)) +
    geom_line(size = line_width) +
    scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +  # Adjust x-axis breaks to show hours
    labs(title = "Smoothed Mean Functions for Pain and No Pain",
         x = "Time (hours)",
         y = "Movement Acceleration (g)",
         color = "Group") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 19, face = "bold"),  # Increase title size
      axis.title = element_text(size = 17),  # Increase axis title size
      axis.text = element_text(size = 15),  # Increase axis text size
      legend.text = element_text(size = 15),  # Increase legend text size
      legend.title = element_text(size = 17)  # Increase legend title size
    )
}


# Function to export the plot and return stats as a data frame row
save_plot_and_return_stats <- function(result, file_name) {
  # Save the plot
  ggsave(filename = file_name, plot = result$plot, width = 4, height = 6, bg = "white")
  # Return the stats as a data frame row
  return(data.frame(
    Comparison = gsub("plot_", "", file_name),
    Group_1 = result$group1_stats$label,
    Group_2 = result$group2_stats$label,
    Group_1_Count = result$group1_stats$count,
    Group_2_Count = result$group2_stats$count,
    Group_1_Mean = result$group1_stats$mean,
    Group_2_Mean = result$group2_stats$mean,
    Group_1_Median = result$group1_stats$median,
    Group_2_Median = result$group2_stats$median
  ))
}


# Function to export the plot and return stats as a data frame row
save_plot_and_return_stats <- function(result, file_name) {
  # Save the plot
  ggsave(filename = file_name, plot = result$plot, width = 4, height = 6, bg = "white")
  # Return the stats as a data frame row
  return(data.frame(
    Comparison = gsub("plot_", "", file_name),
    Group_1 = result$group1_stats$label,
    Group_2 = result$group2_stats$label,
    Group_1_Count = result$group1_stats$count,
    Group_2_Count = result$group2_stats$count,
    Group_1_Mean = result$group1_stats$mean,
    Group_2_Mean = result$group2_stats$mean,
    Group_1_Median = result$group1_stats$median,
    Group_2_Median = result$group2_stats$median
  ))
}

