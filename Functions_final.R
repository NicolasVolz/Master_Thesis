#### packages 

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
library(dplyr)



### create basis functions
create_basis_functions <- function(k = 6, basis_type = "fourier", rangeval = c(0, 5)) {
  if (basis_type == "fourier") {
    basis <- create.fourier.basis(rangeval = rangeval, nbasis = k)
  }
  return(basis)
}

basis = create_basis_functions(k = 6, basis_type = "fourier")

### create beta functions
create_beta_function <- function(t, coeffs, basis) {
  beta_fd <- fd(coeffs, basis)
  beta_values <- eval.fd(t, beta_fd)
  return(beta_values)
}

# Function to check if a matrix is a valid covariance matrix
check_covariance_matrix <- function(cov_matrix) {
  # Check if the matrix is symmetric
  is_symmetric <- all(cov_matrix == t(cov_matrix))
  
  # Calculate the eigenvalues to check for positive semi-definiteness
  eigenvalues <- eigen(cov_matrix)$values
  
  # Check if all eigenvalues are non-negative
  is_positive_semi_definite <- all(eigenvalues >= 0)
  
  # Print the results
  if (is_symmetric && is_positive_semi_definite) {
    print("The matrix is a valid covariance matrix.")
  } else {
    if (!is_symmetric) {
      print("The matrix is not symmetric.")
    }
    if (!is_positive_semi_definite) {
      print("The matrix is not positive semi-definite.")
      print("Eigenvalues:")
      print(eigenvalues)
    }
  }
}



force_positive_semidefinite <- function(matrix) {
  eigen_decomp <- eigen(matrix)
  eigen_decomp$values <- pmax(eigen_decomp$values, 0.00001)  
  return(eigen_decomp$vectors %*% diag(eigen_decomp$values) %*% t(eigen_decomp$vectors))
}
##without force you get eigenvalues like -7.888899e-18

##### create covariance function from covariance matrix
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


####### Function to create random covariance matrix
#DF SMALL-> high correlations, df large -> low correlations
generate_covariance_matrix <- function(dimension, df = 1 , scale_factor = 1) {
  # Step 1: Generate the eigenvalues
  eigenvalues <- rchisq(dimension, df = df) * scale_factor
  # Step 2: random orthogonal Matrix
  orthogonal_matrix <- randortho(dimension)
  
  # Step 3: Construct the covariance matrix
  covariance_matrix <- orthogonal_matrix %*% diag(eigenvalues) %*% t(orthogonal_matrix)
  
  # Ensure symmetry (just in case of any numerical precision issues)
  covariance_matrix <- (covariance_matrix + t(covariance_matrix)) / 2
  
  return(covariance_matrix)
}







#-------------------------------------------------------------------------------
### 1 Simulate Data - one Sample - as Function -> n samples in data frame 
#-------------------------------------------------------------------------------
# Helper function to generate Data
generate_simulated_data <- function(n, n_obs = 100,
                                    beta_0_f, beta_1_f, beta_2_f, 
                                    delta_2,
                                    Treatment, Age, 
                                    error_scale = 1, covariance_structure = NULL) {
  
  # Generate random points in the domain
  t_values <- seq(0, 5, length.out = n_obs)
  
  # generate covariance
  Sigma<- covariance_structure(t_values)
  total_variation <- calculate_trace(Sigma)
  
  # generate true mean 
  beta_0 = beta_0_f(t_values)
  beta_1 = beta_1_f(t_values, delta_2)
  beta_2 = beta_2_f(t_values)
  
  
  # Initialize a list to store the simulated data and components
  sim_data <- list()
  sim_data$t_values <- t_values
  
  for (i in 1:n) {
    # Simulate y_smooth
    y_smooth <- simulate_smooth_component(t_values, Treatment[i], Age[i], beta_0, beta_1, beta_2)
    
    # Simulate y_nonsmooth
    y_nonsmooth <- simulate_nonsmooth_component(t_values)
    
    # Simulate error_function
    error_function <- simulate_error_component(t_values, length_scale, error_scale, Sigma)
    
    # Combine components to get raw_data
    raw_data <- y_smooth + y_nonsmooth + error_function
    
    # Store each component and raw_data in the list
    sim_data[[paste0("Sample_", i)]] <- list(
      raw_data = raw_data,
      y_smooth = y_smooth,
      y_nonsmooth = y_nonsmooth,
      error_function = error_function
    )
  }
  
  return(sim_data)
}

# Helper function to simulate smooth component
simulate_smooth_component <- function(t_values, Treatment, Age, beta_0, beta_1, beta_2) {
  mean_function_smooth <- beta_0 + Treatment * beta_1 + Age * beta_2
  return(mean_function_smooth)
}


# Helper function to simulate nonsmooth component
simulate_nonsmooth_component <- function(t_values, nbasis = 6, basis_type = "fourier", noise_sd = 0.4) {
  # Step 1: Generate the original function as random noise at each time point
  original_function <- rnorm(length(t_values), mean = 0, sd = noise_sd)
  
  # Step 2: Create basis functions for smoothing
  create.fourier.basis(c(min(t_values), max(t_values)), nbasis)
  
  # Step 3: Smooth the original function
  smooth_result <- smooth.basis(t_values, original_function, basis)
  smooth_function <- eval.fd(t_values, smooth_result$fd)
  
  # Step 4: Calculate the nonsmooth component by subtracting the smooth part
  nonsmooth_component <- original_function - smooth_function
  
  return(nonsmooth_component)
}



# Helper function to simulate error component
simulate_error_component <- function(t_values, length_scale, error_scale, covariance_structure = NULL) {
  mean_function_error <- function(t) 0 * t
  Sigma_error <- covariance_structure
  return(error_scale * mvrnorm(n = 1, mu = mean_function_error(t_values), Sigma = Sigma_error))
}


# Calculate the trace of a covariance matrix
calculate_trace <- function(Sigma) {
  return(sum(diag(Sigma)))
}



# Assuming simulated_data is your list
extract_data <- function(sim_data, type = "raw_data") {
  # Valid types
  valid_types <- c("raw_data", "y_smooth", "y_nonsmooth", "error_function")
  
  # Check if the specified type is valid
  if (!(type %in% valid_types)) {
    stop("Invalid type specified. Choose from 'raw_data', 'y_smooth', 'y_nonsmooth', 'error_function'.")
  }
  
  # Extract t_values
  t_values <- sim_data$t_values
  
  # Initialize a data frame to store the desired data type
  data_df <- data.frame(t_values = t_values)
  
  # Loop through each sample in the list
  for (i in seq_along(sim_data)) {
    if (names(sim_data)[i] != "t_values") {
      sample_data <- sim_data[[i]]
      col_name <- paste0("Sample_", i)
      if (type %in% names(sample_data)) {
        data_df[[col_name]] <- sample_data[[type]]
      } else {
        stop(paste("Invalid type specified. Choose from", paste(valid_types, collapse = ", ")))
      }
    }
  }
  
  return(data_df)
}





#-------------------------------------------------------------------------------
### 2. Smooth Data
#-------------------------------------------------------------------------------


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
    rss_values <- numeric(max_basis - 2 + 1)
    
    # Loop over a range of basis numbers to calculate RSS
    for (nbasis in 2:max_basis) {
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
      rss_values[nbasis - 2 + 1] <- rss
    }
    
    # Store the RSS values for the current subject
    rss_values_list[[i - 1]] <- rss_values
    
    # Create and store the Scree Plot for each subject
    nbasis_values <- seq(2, max_basis)
    rss_data <- data.frame(nbasis = nbasis_values, rss = rss_values)
    rss_plot <- ggplot(rss_data, aes(x = nbasis, y = rss)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Scree Plot"),  # Updated to include subject number
           x = "Number of Basis Functions",
           y = "Residual Sum of Squares (RSS)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))  # Centering the title
    
    plot_list[[i - 1]] <- rss_plot  # Store the plot in the list
  }
  
  return(plot_list)  # Return the list of plots
}





smooth_simulated_data <- function(sim_data, nbasis, basis_type = "fourier") {
  # Extract time values
  T <- sim_data$t_values
  
  # Create basis according to the specified type
  if (basis_type == "fourier") {
    basis <- create.fourier.basis(c(0, 5), nbasis)
  } else if (basis_type == "bspline") {
    basis <- create.bspline.basis(c(0, 5), nbasis)
  } else if (basis_type == "constant") {
    basis <- create.constant.basis(c(0, 5), nbasis)
  } else if (basis_type == "monomial") {
    basis <- create.monomial.basis(c(0, 5), nbasis)
  } else {
    stop("Invalid basis_type. Please choose 'fourier', 'bspline', or 'monomial'.")
  }
  
  # Combine all samples into a matrix
  samples_matrix <- as.matrix(sim_data[, -1])
  
  # Perform smoothing on all samples together
  smooth_result <- smooth.basis(T, samples_matrix, basis)
  
  return(smooth_result)
}









################################################################################
############################-- main function---#################################
################################################################################

simulate_2s <- function(n1 = 20, n2 = 20,
                        beta_0, beta_1, beta_2,
                        cov_1, cov_2,
                        n_basis = 7,
                        error_scale = 1,
                        plot = FALSE,
                        basis = FALSE,
                        delta_2 = 0, 
                        Age_1 = Age_1, 
                        Age_2 = Age_2) {
  
  
  optimal_nbasis = 6 ## by lookling at the scree plot
  
  Treatment_1 = c(rep(1,n1))
  Treatment_2 = c(rep(0,n2))

  
  ### generate data 
  data_1 <- generate_simulated_data(n1, n_obs = 100,
                                    beta_0, beta_1, beta_2,
                                    delta_2,
                                    Treatment_1, Age_1,
                                    error_scale = 1, covariance_structure = cov_1)
  data_2 <- generate_simulated_data(n2, n_obs = 100,
                                    beta_0, beta_1, beta_2,
                                    delta_2,
                                    Treatment_2, Age_2,
                                    error_scale = 1, covariance_structure = cov_2)

  
  
  # Combine data for both groups
  simulated_data_1 <- extract_data(data_1, type = "raw_data")
  simulated_data_2 <- extract_data(data_2, type = "raw_data")
  
  # Combine Group 1 and Group 2 data, accounting for unbalanced sizes
  combined_data <- cbind(simulated_data_1, simulated_data_2[,-1])  # Avoid duplicating t_values column
  
  t_values <- simulated_data_1$t_values
  num_samples <- ncol(df) - 1  # Subtract 1 to account for the t_values column
  
  
  ######### plots 
  if (plot){
    # Plot raw data, y_smooth, y_nonsmooth, and error_function for one sample
    png(filename = "plot_sample_data.png", width = 800, height = 600, res = 110)
    plot_1 = plot_sample_data(data_1, sample_num = 1, t_values = t_values)
    print(plot_1)
    dev.off()
    

    # Create an initial plot with the first sample
    png(filename = "plot_multiple_samples.png", width = 800, height = 600, res =110)
    plot_2 = plot_multiple_samples(simulated_data_1, t_values)
    print(plot_2)
    dev.off()
    
  }
  
  #### plot scree plot
  if(basis){
  scree_plots <- determine_optimal_nbasis(combined_data, basis_type = "fourier", max_basis = 20)
  
  # Set the filename and path for the output plot
  output_filename <- "Scree_Simulation.png"
  # Open a PNG device to save the plot
  png(filename = output_filename, width = 800, height = 600, res = 220)
  print(scree_plots[[5]])
  dev.off()
}
  
  
  ## smooth data 
  smoothed_data_1 <- smooth_simulated_data(simulated_data_1, nbasis = optimal_nbasis)
  smoothed_data_1 <- smoothed_data_1$fd
  smoothed_data_2 <- smooth_simulated_data(simulated_data_2, nbasis = optimal_nbasis)
  smoothed_data_2 <- smoothed_data_2$fd
  smoothed_combined_data <- smooth_simulated_data(combined_data, nbasis = optimal_nbasis)
  smoothed_combined_data <- smoothed_combined_data$fd
  
  
  ## plot smoothing example
  if (plot){
    png(filename = "plot_smoothing_examples.png", width = 800, height = 600, res = 110)
    plot_3 = plot_smoothing_example(simulated_data_1[, 2], smoothed_data_1[1], t_values)
    print(plot_3)
    dev.off()
    
    png(filename = "plot_smoothing_examples_points.png", width = 800, height = 600, res = 110)
    plot_smoothing_example(simulated_data_1[, 2], smoothed_data_1[1], t_values, raw_as_points = TRUE)
    dev.off()
  }
  
  
  ## Design_Matrix_X
  n_beta <- 3
  p <- 3
  groupList <- vector("list", n_beta)
  
  sample_headers_1 <- paste0("Sample_", 1:n1)
  sample_headers_2 <- paste0("Sample_", (n1 + 1):(n1 + n2))
  sample_headers <- c(sample_headers_1, sample_headers_2)
  
  # Assign data to each element in the list with different sample sizes
  groupList[[1]] <- setNames(rep(1, n1 + n2), sample_headers)
  groupList[[2]] <- setNames(c(Treatment_1, Treatment_2), sample_headers)
  groupList[[3]] <- setNames(c(Age_1, Age_2), sample_headers)
  groupList[[3]] < scale(groupList[[3]], center = TRUE, scale = FALSE)  # Age
  names(groupList) <- c('Intercept', 'beta_1', 'beta_2')
  

  
  ## set up beta parameter 
  betabasis      = create.fourier.basis(c(0, 5), n_basis)
  betafdPar      = fdPar(betabasis)
  betaList       = vector("list",p)
  for (j in 1:p) betaList[[j]] = betafdPar
  
  
  ## calculate beta
  fRegressList= fRegress(smoothed_combined_data, groupList, betaList)
  
  
  #  extract the estimated regression coefficients and y-values
  betaestList = fRegressList$betaestlist
  Fit   = fRegressList$yhatfd
  groups  = c("Intercept", "beta_1", "beta_2")
  
  
  # Example usage
  if(plot){
    
    ## get mean functions
    mean_function_1 = extract_data(data_1, type = "y_smooth")
    mean_function_2 = extract_data(data_2, type = "y_smooth")
    
    
    ## ?bger verschiedene level von Age mitteln 
    mean_function_1 = rowMeans(mean_function_1[, -1])
    mean_function_2 = rowMeans(mean_function_2[, -1])
    
    
    # Plot the smoothed data and true mean functions
    png(filename = "plot_smoothed_and_true_means.png", width = 800, height = 600, res = 110)
    plot_4 = plot_smoothed_and_true_means(smoothed_data_1, smoothed_data_2, mean_function_1, mean_function_2, t_values)
    print(plot_4)
    dev.off()
    
    
    # Calculate the approximated functions
    approximated_function_1 <- eval.fd(t_values, betaestList[[1]]$fd + betaestList[[2]]$fd + mean((c(Age_1, Age_2))) * betaestList[[3]]$fd)
    approximated_function_2 <- eval.fd(t_values, betaestList[[1]]$fd + mean((c(Age_1, Age_2))) * betaestList[[3]]$fd)
  }
  
  
  
  
  

  
  ############################# Test Konrad ####################################
  
  # get XX^-1
  X = matrix(c(c(rep(1,(n1+n2))), c(Treatment_1, Treatment_2), c(Age_1, Age_2)),
             byrow = FALSE, ncol = 3)
  XtX_inv <- solve(t(X) %*% X)
  XtX_inv_jj <- XtX_inv[2, 2]  
  
  ## get beta_1
  beta_fd <- betaestList[[2]]$fd
  
  ## get estimated covariance Matrix gamma_hat 
  smoothed_combined_hat = fRegressList$yhatfdobj
  residuals_hat = smoothed_combined_data - smoothed_combined_hat
  gamma_hat = var.fd(residuals_hat)
  
  
  ## calcualte V_hat
  bj_hat = c(beta_fd$coefs)
  residuals_matrix = t(as.matrix(residuals_hat$coefs))
  V_hat = cov(residuals_matrix)
  Tn_K =c(c(bj_hat) %*% solve(V_hat)/(XtX_inv_jj)) %*% c(bj_hat)
  significant_K = Tn_K > qchisq(0.95, df = n_basis)    
  
  
  
  ## with V_hat compute estimated covariance function: 
  covariance_structure_estimated = function(t_values){
    create_covariance_structure(t_values, V_hat)
  }
  Sigma_estimated = covariance_structure_estimated(t_values)
  

  ## true mean and covariance converge towards values given?
  if (plot){
    
    ##b_values
    beta_fd_0 <- betaestList[[1]]$fd 
    beta_fd_1 <- betaestList[[2]]$fd  
    beta_fd_2 <- betaestList[[3]]$fd
    
    bj_hat_0 = c(beta_fd_0$coefs)
    bj_hat_1 = c(beta_fd_1$coefs)
    bj_hat_2 = c(beta_fd_2$coefs)
    
    XtX_inv_diag <- diag(XtX_inv)
    
    # Compute standard errors
    se_bj_hat_0 <- sqrt(diag(V_hat * XtX_inv_diag[1]))
    se_bj_hat_1 <- sqrt(diag(V_hat * XtX_inv_diag[2]))
    se_bj_hat_2 <- sqrt(diag(V_hat * XtX_inv_diag[3]))
  
    
    # Construct 95% confidence intervals
    ci_bj_hat_0 <- cbind(bj_hat_0 - 1.96 * se_bj_hat_0, bj_hat_0 + 1.96 * se_bj_hat_0)
    ci_bj_hat_1 <- cbind(bj_hat_1 - 1.96 * se_bj_hat_1, bj_hat_1 + 1.96 * se_bj_hat_1)
    ci_bj_hat_2 <- cbind(bj_hat_2 - 1.96 * se_bj_hat_2, bj_hat_2 + 1.96 * se_bj_hat_2)
    
    
    png(filename = "visualize_all_vector_differences.png", width = 800, height = 600, res = 110)
    plot_5 = visualize_all_vector_differences(bj_hat_0, beta_0_v, bj_hat_1, 
                                     beta_1_v, bj_hat_2, beta_2_v, 
                                     ci_bj_hat_0, ci_bj_hat_1, ci_bj_hat_2)
    print(plot_5)
    dev.off()
    
    
    # V matrices
    
    rownames(Sigma_k1) <- rownames(V_hat)
    colnames(Sigma_k1) <- colnames(V_hat)
    
    VV = diag(diag(Sigma_k1))
    rownames(VV) <- rownames(V_hat)
    colnames(VV) <- colnames(V_hat)

    
    png(filename = "plot_V_estimated_1.png", width = 800, height = 600, res = 110)
    plot_7 = plot_V_matrix(V_hat, title = expression(~bold(hat(V))))
    print(plot_7)
    dev.off()
    
    print(round(V_hat,2))
    
    ## covariance
    Sigma_true_1 <- cov_1(t_values)
    Sigma_true_2 <- cov_2(t_values)
    
    #correlation
    
    cor_1_true = cov2cor(Sigma_true_1)
    cor_2_true = cov2cor(Sigma_true_1)
    cor_es = cov2cor(Sigma_estimated)
    
    
    
    png(filename = "plot_covariance_estimated.png", width = 800, height = 600, res = 110)
    plot_10 = plot_covariance_matrix(Sigma_estimated, bquote(hat(gamma)("t", "t'")))
    print(plot_10)
    dev.off()
    
    
    plot_covariance_matrix(cor_1_true, bquote(cor("t", "t'")))
    plot_covariance_matrix(cor_2_true, bquote(cor("t", "t'")))
    plot_covariance_matrix(cor_es, bquote(widehat(cor)("t", "t'")))
    
    
    
    beta_0_hat = as.numeric(eval.fd(betaestList[[1]]$fd, t_values))
    beta_1_hat = as.numeric(eval.fd(betaestList[[2]]$fd, t_values))
    beta_2_hat = as.numeric(eval.fd(betaestList[[3]]$fd, t_values))
    
    
    ### Apply the function to each beta function
    se_beta_0_hat <- pointwise_sd(beta_fd_0, XtX_inv_diag[1], gamma_hat, t_values)
    se_beta_1_hat <- pointwise_sd(beta_fd_1, XtX_inv_diag[2], gamma_hat, t_values)
    se_beta_2_hat <- pointwise_sd(beta_fd_2, XtX_inv_diag[3], gamma_hat, t_values)
    
    # Construct pointwise confidence intervals
    ci_bj_hat_0 <- cbind(beta_0_hat - 1.96 * se_beta_0_hat, beta_0_hat + 1.96 * se_beta_0_hat)
    ci_bj_hat_1 <- cbind(beta_1_hat - 1.96 * se_beta_1_hat, beta_1_hat + 1.96 * se_beta_1_hat)
    ci_bj_hat_2 <- cbind(beta_2_hat - 1.96 * se_beta_2_hat, beta_2_hat + 1.96 * se_beta_2_hat)
  
    
    # Plot beta with confidence intervals
    
    png(filename = "plot_true_vs_approximated_1.png", width = 800, height = 600, res = 110)
    plot_11 = plot_true_vs_approximated(t_values, 
                              beta_0_f(t_values), 
                              as.numeric(eval.fd(betaestList[[1]]$fd, t_values)),
                              approx_color = "blue",
                              title = expression(beta[0](t)), 
                              true_label = expression(beta[0](t)),
                              approx_label = expression(hat(beta)[0](t)),
                              ci_lower = ci_bj_hat_0[, 1],  # Lower bound of CI
                              ci_upper = ci_bj_hat_0[, 2]   # Upper bound of CI
    )
    print(plot_11)
    dev.off()
    
    png(filename = "plot_true_vs_approximated_2.png", width = 800, height = 600, res = 110)
    plot_12 = plot_true_vs_approximated(t_values, 
                              beta_1_f(t_values, delta_2), 
                              as.numeric(eval.fd(betaestList[[2]]$fd, t_values)),
                              approx_color = "blue",
                              title = expression(beta[1](t)), 
                              true_label = expression(beta[1](t)),
                              approx_label = expression(hat(beta)[1](t)),
                              ci_lower = ci_bj_hat_1[, 1],  # Lower bound of CI
                              ci_upper = ci_bj_hat_1[, 2]   # Upper bound of CI
    )
    print(plot_12)
    dev.off()
    
    
    png(filename = "plot_true_vs_approximated_3.png", width = 800, height = 600, res = 110)
    plot_13 = plot_true_vs_approximated(t_values, 
                              beta_2_f(t_values), 
                              as.numeric(eval.fd(betaestList[[3]]$fd, t_values)),
                              approx_color = "blue",
                              title = expression(beta[2](t)), 
                              true_label = expression(beta[2](t)),
                              approx_label = expression(hat(beta)[2](t)),
                              ci_lower = ci_bj_hat_2[, 1],  # Lower bound of CI
                              ci_upper = ci_bj_hat_2[, 2]   # Upper bound of CI
    )
    print(plot_13)
    dev.off()

  }
  
  
  
  #############################------L2 Test---------###########################

  # Calculate the covariance matrix for b_jh
  Cov_bjh_bjh <- XtX_inv_jj * V_hat
  
  # Perform SVD
  svd_results <- svd(Cov_bjh_bjh)
  
  # Extract eigenvalues (singular values squared give the eigenvalues for covariance matrices)
  eigenvalues <- svd_results$d
  
  # chi squared mixture
  chi_squares <- sapply(eigenvalues, function(lambda) {
    rchisq(1, df = 1) * lambda
  })
  
  
  ### shatterwhite welch approximation
  gamma = Cov_bjh_bjh
  
  # Compute gamma^otimes2 - you may need to define this integral computation
  gamma_otimes2 <- compute_gamma_otimes2(gamma)  # Assume this function computes the integral and returns a matrix
  
  # Calculate traces
  trace_gamma <- sum(diag(gamma))
  trace_gamma_otimes2 <- sum(diag(gamma_otimes2))
  
  # Calculate beta and k using Welch-Satterthwaite
  beta_hat_k <- trace_gamma_otimes2 / trace_gamma
  k_hat_k <- (trace_gamma^2) / trace_gamma_otimes2
  
  # Usage in statistical tests (if you need to compute critical values or similar)
  critical_value_2 <- beta_hat_k * qchisq(0.95, df = k_hat_k)
  
  
  # Tn_K from actual data might look like this:
  # Assuming XtX_inv_jj is calculated and V_hat is derived from actual data
  Tn_K_actual <- (t(bj_hat) %*% bj_hat)
  
  # Test for significance
  #significant_K_actual <- Tn_K_actual > critical_value
  significant_K_L2 <- Tn_K_actual > critical_value_2
  
  
  ##### plot approximation through simulation
  if (plot){
    
    ### direct approximation
    n_sims <- 10000
    # Simulate the test statistic under the null hypothesis
    simulated_stats <- replicate(n_sims, sum(sapply(eigenvalues, function(lambda) rchisq(1, df = 1) * lambda)))
    critical_value <- quantile(simulated_stats, 0.95)
    ## plot approximated distributions
    # Generate data for Welch-Satterthwaite distribution
    ws_distribution <- rchisq(n_sims, df = k_hat_k) * beta_hat_k
    
    # Plot both distributions for comparison
    hist(simulated_stats, breaks = 50, probability = TRUE, main = "Comparison of Approximated Distributions L2",
         xlab = "Test Statistic Value", col = rgb(0.2, 0.4, 0.6, 0.7), border = "white")
    lines(density(ws_distribution), col = "red", lwd = 2)
    
    legend("topright", legend = c("Simulated Distribution", "Welch-Satterthwaite Approximation"),
           col = c(rgb(0.2, 0.4, 0.6, 0.7), "red"), lwd = 2, bty = "n")
    abline(v = critical_value, col = "blue", lwd = 2, lty = 2)
    abline(v = critical_value_2, col = "green", lwd = 2, lty = 2)
    legend("right", legend = c("Critical Value (Simulated)", "Critical Value (WS Approx)"),
           col = c("blue", "green"), lwd = 2, lty = 2, bty = "n")
    
  }
  
  
  #############################-------F-Test-------- ###########################
  ## nenner
  nenner = Tn_K_actual
  
  # Zahler
  zahler = trace_gamma
  
  ## test statistic
  Fn = (nenner) / (zahler)
  
  ## approximate Distribution
  df1_k = k_hat_k
  df2_k = ((n1+n2)-p-1)*k_hat_k
  
  critical_value_3 = qf(0.95, df1 = df1_k, df2 = df2_k)
  significant_F_k = Fn > critical_value_3
  
  
  
  ##### plot approximation through simulation
  if (plot){
    
    ### direct approximation
    n_sims <- 10000
    # Simulate the test statistic under the null hypothesis
    simulated_stats_1 <- replicate(n_sims, sum(sapply(eigenvalues, function(lambda) rchisq(1, df = 1) * lambda)))
    simulated_stats_2 <- replicate(n_sims, sum(sapply(eigenvalues, function(lambda) rchisq(1, df = (n1+n2)-p-1) * lambda)))/(n1+n2-p-1)
    simulated_stats_F = simulated_stats_1 / simulated_stats_2
    critical_value <- quantile(simulated_stats_F, 0.95)
    
    
    ## plot approximated distributions
    # Generate data for Welch-Satterthwaite distribution
    ws_distribution <- rf(n_sims,df1 = df1_k, df2 = df2_k)
    # Plot both distributions for comparison
    hist(simulated_stats_F, breaks = 50, probability = TRUE, main = "Comparison of Approximated Distributions F",
         xlab = "Test Statistic Value", col = rgb(0.2, 0.4, 0.6, 0.7), border = "white")
    lines(density(ws_distribution), col = "red", lwd = 2)
    
    legend("topright", legend = c("Simulated Distribution", "Welch-Satterthwaite Approximation"),
           col = c(rgb(0.2, 0.4, 0.6, 0.7), "red"), lwd = 2, bty = "n")
    abline(v = critical_value, col = "blue", lwd = 2, lty = 2)
    abline(v = critical_value_3, col = "green", lwd = 2, lty = 2)
    legend("right", legend = c("Critical Value (Simulated)", "Critical Value (WS Approx)"),
           col = c("blue", "green"), lwd = 2, lty = 2, bty = "n")
    
  }
  
  
  
  
  
  ##########################----Pointwise Test------- ##########################
  
  if(plot){
    
    ## fixed_t
    t = 4
    points_group_1 = eval.fd(smoothed_data_1, t)
    points_group_2 = eval.fd(smoothed_data_2, t)
    #print(t.test(points_group_1, points_group_2, var.equal = TRUE)$conf.int)
    
    ## confidence interval
    # Create t values
    t_values <- seq(0, 5, length.out = 100)
    
    # Initialize vectors to store confidence interval bounds
    lower_bounds <- numeric(length(t_values))
    upper_bounds <- numeric(length(t_values))
    
    # Loop through each t value
    for (i in seq_along(t_values)) {
      t <- t_values[i]
      
      # Perform t-test for each t value
      points_group_1 <- eval.fd(smoothed_data_1, t)
      points_group_2 <- eval.fd(smoothed_data_2, t)
      t_test_result <- t.test(points_group_1, points_group_2, var.equal = TRUE)
      
      # Extract confidence interval bounds
      lower_bound <- t_test_result$conf.int[1]
      upper_bound <- t_test_result$conf.int[2]
      
      # Store confidence interval bounds
      lower_bounds[i] <- lower_bound
      upper_bounds[i] <- upper_bound
    }
    
    # Plot beta_2(t) with confidence band
    plot(t_values, beta_1(t_values, delta_2), type = "l", col = "black", ylim = c(min(lower_bounds), max(upper_bounds)),
         main = "Beta Function with Confidence Band", xlab = "Time", ylab = "Value")
    polygon(c(t_values, rev(t_values)), c(lower_bounds, rev(upper_bounds)), col = "lightblue", border = NA)
    lines(t_values, beta_1(t_values, delta_2), type = "l", col = "black")
    lines(betaestList[[2]]$fd, col = "blue")
    abline(h = 0, col = "red", lty = 2)
    legend("topright", legend = c("True Beta Function", "Confidence Band", "Estimated Beta Function"), 
           col = c("black", "lightblue", "blue"), lty = c(1, 1, 1))
    
    ## permutation Tests 
    tperm.fd(smoothed_data_1, smoothed_data_2)
    # Fperm.fd(smoothed_combined_data, groupList, betaList)
  }
  
  
  

  return (c(significant_K, significant_F_k, significant_K_L2))
  
}



###############----------Helper Functiions----##################################

# Calculate pointwise standard deviations for all beta_j
pointwise_sd <- function(beta_fd, XtX_inv_diag, gamma_hat, t_values) {
  
  # Calculate the pointwise variances
  pointwise_variances <- sapply(1:length(t_values), function(i) {
    t_val <- t_values[i]
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



plot_simulation_results <- function(df_melted, title, delta_2, include_legend = FALSE, common_y_limit = c(0, 1)) {
  # Handle both n1 and n2 for plotting
  df_melted$n_combined <- df_melted$n1 + df_melted$n2  # Combine n1 and n2 if needed
  
  max_n <- max(df_melted$n_combined, na.rm = TRUE)  # Ensure max_n is a finite number
  
  if (!is.finite(max_n) || max_n <= 0) {
    stop("Error: Unable to determine a valid range for 'n_combined' in the plot.")
  }
  
  # Conditionally determine if y-axis labels should be shown
  show_y_axis_labels <- delta_2 %in% c(0, 1)
  
  p <- ggplot(df_melted, aes(x = n_combined, y = Power, color = Test, shape = Test)) +
    geom_line(aes(linetype = Test), size = 1) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = seq(0, max_n, by = 50)) +
    scale_y_continuous(limits = common_y_limit, breaks = seq(0, 1, by = 0.1)) +  # Common y-axis limit
    labs(
      title = title,
      x = NULL,  # Remove x-axis label
      y = NULL   # Remove y-axis label
    ) +
    theme_minimal() +
    scale_color_manual(values = c("L2 Test" = "blue", "F Test" = "red", "K Test" = "orange")) +
    scale_shape_manual(values = c("L2 Test" = 18, "F Test" = 17, "K Test" = 16)) +
    scale_linetype_manual(values = c("L2 Test" = "dotted", "F Test" = "dashed", "K Test" = "solid")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22),  # Increase title size
      axis.text.x = element_text(size = 16),  # Increase x-axis text size
      axis.text.y = if (show_y_axis_labels) element_text(size = 16) else element_blank(),  # Conditionally show y-axis text
      axis.title.x = element_text(size = 16),  # Increase x-axis title size
      axis.title.y = if (show_y_axis_labels) element_text(size = 16) else element_blank(),  # Conditionally show y-axis title
      legend.title = element_blank(),
      legend.position = if (include_legend) "right" else "none",
      legend.text = element_text(size = 18)  # Increase legend text size
    )
  
  if (delta_2 == 0) {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dashed", color = "black")
  }
  
  return(p)
}



create_decreasing_vector <- function(size, initial_value = 1, noise_sd = 1) {
  values <- sapply(1:size, function(k) initial_value * (1/(k^2)) + rnorm(1, mean = 0, sd = noise_sd))
  return(values)
}


# Function to extract the evaluation data from fd objects
extract_fd_data <- function(fd_object, t_values) {
  eval_fd <- eval.fd(t_values, fd_object)
  data.frame(t = rep(t_values, ncol(eval_fd)),
             value = as.vector(eval_fd),
             curve = factor(rep(1:ncol(eval_fd), each = length(t_values))))
}








################################################################################
####################################------PLOTS---------######################## 
################################################################################



### general theme
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
  legend.position = "right",  # Position legend to the right
  strip.text = element_text(face = "bold", size = 22),  # Facet strip text size
  strip.background = element_blank(),  # Remove background box
  panel.spacing = unit(0.5, "lines"),  # Space between panels
  panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Single line border around each panel
)




visualize_all_vector_differences <- function(bj_hat_0, beta_0_v, bj_hat_1, beta_1_v, bj_hat_2, beta_2_v, ci_bj_hat_0, ci_bj_hat_1, ci_bj_hat_2, title = "Comparison of Estimated and True Coefficients") {
  
  bj_hat_2 = bj_hat_2 * 50
  beta_2_v = beta_2_v * 50
  ci_bj_hat_2 = ci_bj_hat_2 * 50
  
  # Ensure matching lengths between coefficients and confidence intervals
  if (length(bj_hat_0) != nrow(ci_bj_hat_0) || length(bj_hat_1) != nrow(ci_bj_hat_1) || length(bj_hat_2) != nrow(ci_bj_hat_2)) {
    stop("Mismatch in lengths between coefficients and confidence intervals.")
  }
  
  # Create a combined data frame for ggplot
  data <- data.frame(
    Index = rep(seq(0, length(bj_hat_0) - 1), times = 3),  # Start Index from 0
    Beta_Hat = c(bj_hat_0, bj_hat_1, bj_hat_2),
    Beta_V = c(beta_0_v, beta_1_v, beta_2_v),
    Lower_CI = c(ci_bj_hat_0[, 1], ci_bj_hat_1[, 1], ci_bj_hat_2[, 1]),
    Upper_CI = c(ci_bj_hat_0[, 2], ci_bj_hat_1[, 2], ci_bj_hat_2[, 2]),
    Coefficient = factor(rep(c("beta_0", "beta_1", "beta_2"), each = length(bj_hat_0)))
  )
  
  # Convert Coefficient labels to LaTeX format using expression()
  data$Coefficient <- factor(data$Coefficient,
                             levels = c("beta_0", "beta_1", "beta_2"),
                             labels = c(expression(bold(b[0])), expression(bold(b[1])), expression(bold(b[2]))))
  
  # Create the plot
  p <- ggplot(data, aes(x = Index)) +
    geom_point(aes(y = Beta_Hat, color = "Estimated~(hat(b))"), size = 2, shape = 16) +  # Estimated points
    geom_point(aes(y = Beta_V, color = "True~(b)"), size = 2.5, shape = 17) +  # True points
    geom_linerange(aes(ymin = Lower_CI, ymax = Upper_CI), color = "red") +  # Add CI error bars
    facet_wrap(~ Coefficient, scales = "fixed", labeller = label_parsed) +
    labs(
      title = title,
      x = expression(Index~h),  # X-axis label with "Index h"
      y = "Value"
    ) +
    scale_x_continuous(breaks = 0:(length(bj_hat_0) - 1), limits = c(0, length(bj_hat_0) - 1)) +  # Set x-axis to start from 0
    scale_color_manual(
      name = "",  # Remove the legend title
      values = c("Estimated~(hat(b))" = "red", "True~(b)" = "green"),
      labels = c(expression(hat(b)[h]), expression(b[h]))  # Correctly display in legend
    ) +
    theme_minimal() +
    custom_theme + 
    theme(panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_line(size = 0.25))
  
  # Print the plot
  print(p)
  return(p)
}




# Function to plot raw data, y_smooth, y_nonsmooth, and error_function for one sample
plot_sample_data <- function(data, sample_num = 1, t_values) {
  sample_data <- data[[paste0("Sample_", sample_num)]]
  raw_data <- sample_data$raw_data
  y_smooth <- sample_data$y_smooth
  y_nonsmooth <- sample_data$y_nonsmooth
  error_function <- sample_data$error_function
  
  # Convert data to a data frame for ggplot
  plot_data <- data.frame(
    t = rep(t_values, 4),
    value = c(raw_data, y_smooth, y_nonsmooth, error_function),
    type = factor(rep(c("y_raw", "y_smooth", "y_nonsmooth", "error_function"), each = length(t_values)),
                  levels = c("y_raw", "y_smooth", "y_nonsmooth", "error_function"))
  )
  
  # Plot the data using ggplot2 for consistency
  plot <- ggplot(plot_data, aes(x = t, y = value, color = type, linetype = type)) +
    geom_line(size = 1) +
    scale_color_manual(
      values = c("y_raw" = "red", "y_smooth" = "blue", "y_nonsmooth" = "green", "error_function" = "purple"),
      labels = c(expression(f[i](t)), expression(bold(X)[i*"."  ] * bold(beta)(t)), expression(y[i]^(ns)~(t)), expression(epsilon[i](t)))
    ) +
    scale_linetype_manual(
      values = c("y_raw" = "solid", "y_smooth" = "solid", "y_nonsmooth" = "solid", "error_function" = "solid"),
      labels = c(expression(f[i](t)), expression(bold(X)[i*"." ] * bold(beta)(t)), expression(y[i]^(ns)~(t)), expression(epsilon[i](t)))
    ) +
    labs(title = "Example", x = "t", y = "Value", color = "Legend", linetype = "Legend") +
    theme_minimal() +
    scale_x_continuous(
      breaks = seq(min(t_values), max(t_values), length.out = 13),  # Adjust based on the data range
      labels = seq(0, 24, by = 2)  # Display labels from 0 to 24 at intervals of 2
    )+
  
    custom_theme
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}




# Plotting function for covariance matrices
plot_covariance_matrix <- function(Sigma, title, x = "t" , y = "t'") {
  Sigma_melted <- melt(Sigma)
  plot <- ggplot(Sigma_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                         limit = c(min(Sigma_melted$value), max(Sigma_melted$value))) +
    scale_x_continuous(
      breaks = seq(min(Sigma_melted$Var1), max(Sigma_melted$Var1), length.out = 13), # 13 evenly spaced ticks
      labels = seq(0, 24, by = 2)  # Labels at intervals of 2 (0, 2, 4, ..., 24)
    ) +
    scale_y_continuous(
      breaks = seq(min(Sigma_melted$Var2), max(Sigma_melted$Var2), length.out = 13), # 13 evenly spaced ticks
      labels = seq(0, 24, by = 2)  # Labels at intervals of 2 (0, 2, 4, ..., 24)
    ) +
    labs(title = title, x = x, y = y) +
    custom_theme +
    theme(
      panel.border = element_blank(),  # Remove the panel border
      axis.line = element_blank()  # Remove axis lines to avoid the box
    )
  
  print(plot)  # Ensure the plot is printed
  return(plot)
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
    custom_theme +
    theme(
      panel.border = element_blank(),  # Remove the panel border
      axis.line = element_blank()  # Remove axis lines to avoid the box
    )
  if (!is_numeric) {
    plot <- plot +
      scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5", "6")) +  # Change X-axis labels
      scale_y_discrete(labels = c("0", "1", "2", "3", "4", "5", "6"))    # Change Y-axis labels
  }
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}



# Function to plot smoothed data and true mean functions using ggplot2
plot_smoothed_and_true_means <- function(smoothed_data_1, smoothed_data_2, mean_function_1, mean_function_2, t_values) {
  # Extract data for plotting
  smoothed_data_1_eval <- extract_fd_data(smoothed_data_1, t_values)
  smoothed_data_2_eval <- extract_fd_data(smoothed_data_2, t_values)
  
  # Prepare data for ggplot2
  smoothed_data_1_eval$type <- "smoothed_data_1"
  smoothed_data_2_eval$type <- "smoothed_data_2"
  
  mean_data <- data.frame(
    t = rep(t_values, 2),
    value = c(mean_function_1, mean_function_2),
    curve = factor(rep(c("true_mean_1", "true_mean_2"), each = length(t_values))),
    type = factor(rep(c("true_mean_1", "true_mean_2"), each = length(t_values)))
  )
  
  # Combine all data
  plot_data <- bind_rows(smoothed_data_1_eval, smoothed_data_2_eval, mean_data)
  
  # Plot using ggplot2
  plot <- ggplot(plot_data, aes(x = t, y = value, group = interaction(curve, type), color = type, linetype = type, size = type)) +
    geom_line() +
    scale_color_manual(
      values = c(
        "smoothed_data_1" = "red", 
        "smoothed_data_2" = "black", 
        "true_mean_1" = "red", 
        "true_mean_2" = "black"
      ),
      labels = c(
        expression(y[i1](t)), 
        expression(y[i2](t)), 
        expression(mu[1](t)), 
        expression(mu[2](t))
      )
    ) +
    scale_linetype_manual(
      values = c(
        "smoothed_data_1" = "dashed", 
        "smoothed_data_2" = "dashed", 
        "true_mean_1" = "solid", 
        "true_mean_2" = "solid"
      ),
      labels = c(
        expression(y[i1](t)), 
        expression(y[i2](t)), 
        expression(mu[1](t)), 
        expression(mu[2](t))
      )
    ) +
    scale_size_manual(
      values = c(
        "smoothed_data_1" = 0.5, 
        "smoothed_data_2" = 0.5, 
        "true_mean_1" = 2, 
        "true_mean_2" = 2
      ),
      labels = c(
        expression(y[i1](t)), 
        expression(y[i2](t)), 
        expression(mu[1](t)), 
        expression(mu[2](t))
      )
    ) +
    labs(
      title = "Smoothed data and true mean functions", 
      x = "t", 
      y = "Value", 
      color = "Legend", 
      linetype = "Legend", 
      size = "Legend"
    ) +
    scale_x_continuous(
      breaks = seq(min(t_values), max(t_values), length.out = 13),  # Adjust based on the data range
      labels = seq(0, 24, by = 2)  # Display labels from 0 to 24 at intervals of 2
    )+
    custom_theme
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}





plot_smoothing_example <- function(simulated_data, smoothed_data, t_values, raw_as_points = FALSE) {
  # Convert data to a long format data frame
  plot_data <- data.frame(
    t = rep(t_values, 2),
    value = c(simulated_data, eval.fd(smoothed_data, t_values)),
    type = factor(rep(c("raw_data", "smoothed_data"), each = length(t_values)),
                  levels = c("raw_data", "smoothed_data"))
  )
  
  # Start the ggplot
  plot <- ggplot(plot_data, aes(x = t, y = value, color = type))
  
  # Add raw data as points or line based on raw_as_points argument
  if (raw_as_points) {
    plot <- plot + geom_point(data = subset(plot_data, type == "raw_data"), size = 2)  # Plot raw data as points
  } else {
    plot <- plot + geom_line(data = subset(plot_data, type == "raw_data"), size = 1)  # Plot raw data as line
  }
  
  # Add the smoothed data as a line
  plot <- plot + geom_line(data = subset(plot_data, type == "smoothed_data"), size = 1) +
    scale_color_manual(
      values = c("raw_data" = "black", "smoothed_data" = "red"),
      labels = c(expression(f[i](t)), expression(y[i]~(t)))
    ) +
    # Handle linetype depending on raw_as_points
    scale_linetype_manual(
      values = c("raw_data" = ifelse(raw_as_points, "blank", "solid"), "smoothed_data" = "solid"),  
      guide = 'none'  # Disable linetype legend
    ) +
    scale_x_continuous(
      breaks = seq(min(t_values), max(t_values), length.out = 13),  # Adjust based on the data range
      labels = seq(0, 24, by = 2)  # Display labels from 0 to 24 at intervals of 2
    ) +    
    
    
    guides(color = guide_legend(override.aes = list(
      linetype = c(ifelse(raw_as_points, "blank", "solid"), "solid"), 
      shape = c(ifelse(raw_as_points, 16, NA), NA)))) +  # Shape: point for raw data if raw_as_points = TRUE
    labs(title = "Smoothing Example", x = "t", y = "Value", color = "Legend") +
    custom_theme
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}






# Function to plot true mean functions and estimates with LaTeX labels and confidence intervals
plot_true_vs_approximated <- function(t_values, mean_function, approximated_function, approx_color, title, true_label, approx_label, ci_lower = NULL, ci_upper = NULL) {
  
  # Prepare data for ggplot2
  plot_data <- data.frame(
    t = rep(t_values, 2),
    value = c(mean_function, approximated_function),
    type = factor(rep(c("True", "Approximated"), each = length(t_values)))
  )
  
  # Create a named vector for labels to ensure the correct mapping
  label_map <- c("True" = true_label, "Approximated" = approx_label)
  
  # Create the base plot
  plot <- ggplot(plot_data, aes(x = t, y = value, color = type, linetype = type, size = type)) +
    geom_line() +
    scale_color_manual(values = c("True" = "black", "Approximated" = approx_color),
                       labels = label_map) +
    scale_linetype_manual(values = c("True" = "solid", "Approximated" = "dashed"),
                          labels = label_map) +
    scale_size_manual(values = c("True" = 1, "Approximated" = 1),
                      labels = label_map) +
    labs(title = title, x = "t", y = "Value", color = "Legend", linetype = "Legend", size = "Legend") +
    custom_theme +
    geom_line(aes(x = t, y = value), data = subset(plot_data, type == "True"), size = 1, linetype = "solid") +
    geom_line(aes(x = t, y = value), data = subset(plot_data, type == "Approximated"), size = 1, linetype = "dashed") +
    theme(legend.title = element_blank())+  # Remove the legend title
  
  scale_x_continuous(
    breaks = seq(min(t_values), max(t_values), length.out = 13),  # Adjust based on the data range
    labels = seq(0, 24, by = 2)  # Display labels from 0 to 24 at intervals of 2
  )
  
  # Add confidence intervals if provided
  if (!is.null(ci_lower) & !is.null(ci_upper)) {
    ci_data <- data.frame(
      t = t_values,
      lower = ci_lower,
      upper = ci_upper
    )
    
    plot <- plot +
      geom_ribbon(aes(x = t, ymin = lower, ymax = upper), data = ci_data, fill = approx_color, alpha = 0.2, inherit.aes = FALSE)
  }
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}








# Function to plot multiple samples with ggplot2
plot_multiple_samples <- function(simulated_data, t_values) {
  # Convert simulated_data to a long format data frame
  plot_data <- data.frame(
    t = rep(t_values, ncol(simulated_data) - 1),
    value = unlist(simulated_data[,-1]),
    sample = factor(rep(2:ncol(simulated_data), each = length(t_values)))
  )
  
  # Plot the data using ggplot2 for consistency
  plot <- ggplot(plot_data, aes(x = t, y = value, group = sample, color = sample)) +
    geom_line(size = 1) +
    scale_color_manual(values = sample(colors(), ncol(simulated_data) - 1, replace = TRUE)) +
    labs(title = "All Samples", x = "t", y = "Value") +
    scale_x_continuous(
      breaks = seq(min(t_values), max(t_values), length.out = 13),  # Adjust based on the data range
      labels = seq(0, 24, by = 2)  # Display labels from 0 to 24 at intervals of 2
    )+
    custom_theme + 
    theme(
      legend.position = "none"  # No legend
    )
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}







# Plotting function for V matrices with adjusted axis labels
plot_V_matrix_alternative <- function(V_matrix, title, x = expression(epsilon[h]), y = expression(epsilon[h])) {
  library(ggplot2)
  library(reshape2)
  
  V_melted <- melt(V_matrix)
  
  # Adjust scale accordingly
  plot <- ggplot(V_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                         limit = c(min(V_melted$value), max(V_melted$value))) +
    labs(title = title, x = x, y = y) +
    custom_theme +
    theme(
      panel.border = element_blank(),  # Remove the panel border
      axis.line = element_blank()      # Remove axis lines to avoid the box
    ) +
    scale_x_continuous(breaks = 1:7, labels = as.character(0:6)) +  # Set X-axis labels from 0 to 6
    scale_y_continuous(breaks = 1:7, labels = as.character(0:6))    # Set Y-axis labels from 0 to 6
  
  print(plot)  # Ensure the plot is printed
  return(plot)
}






# Function to create individual plots with labels in LaTeX
plot_single_basis <- function(x, y, label) {
  data <- data.frame(x = x, y = y)
  ggplot(data, aes(x = x, y = y)) +
    geom_line(size = 1) +
    theme_minimal(base_size = 14) +
    labs(title = bquote(.(label)), x = "t", y = "") +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = -1, size = 12, face = "bold"), # Center and style the label
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5, 5, 5, 5),
      panel.background = element_rect(fill = "gray95", color = NA), # Set background color
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}












