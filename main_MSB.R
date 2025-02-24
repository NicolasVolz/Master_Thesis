#### Main Run ####


### define pain ids

# pain_ids_chron
# pain_ids_kroff_40
# pain_ids_kroff_50
# pain_ids_beeintraechtigung

pain_ids = pain_ids_chron
no_pain_ids = no_pain_ids_chron



####################### one Data Set Example ###################################

setwd("/Users/nicol/Desktop/SP5/data/ID Results für Nicolas Volz")
file_path <- "Results.csv"

# Step 1: Read data
data_wide <- read_and_transform_data(file_path, xlsx = FALSE)

scree_plots <- determine_optimal_nbasis(data_wide)


setwd("/Users/nicol/Desktop/SP5/Results/")
# To display a specific plot:
print(scree_plots[[5]])  # Print the second plot, etc.


# Set the filename and path for the output plot
output_filename <- "Scree.png"
# Open a PNG device to save the plot
png(filename = output_filename, width = 800, height = 800, res = 180)
print(scree_plots[[5]])
dev.off()

##optimal 25 basis functions


# Step 2: Smooth data and plot example
smoothed_data_orig_2 <- smooth_data(data_wide, nbasis = 10, basis_type = "fourier")
smoothed_data_orig <- smooth_data(data_wide, nbasis = 25, basis_type = "fourier")
smoothed_data_orig_3 <- smooth_data(data_wide, nbasis = 100, basis_type = "fourier")


# Set the filename and path for the output plot
output_filename <- "raw_and_smoothed_data_day1_10_Basis.png"
# Open a PNG device to save the plot
png(filename = output_filename, width = 800, height = 800, res = 140)
plot_raw_and_smoothed_data(data_wide$MovementAcceleration_g_d5, smoothed_data_orig_2, day = 5)
dev.off()

# Set the filename and path for the output plot
output_filename <- "raw_and_smoothed_data_day1_25_Basis.png"
# Open a PNG device to save the plot
png(filename = output_filename, width = 800, height = 800, res = 140)
plot_raw_and_smoothed_data(data_wide$MovementAcceleration_g_d5, smoothed_data_orig, day = 5)
dev.off()

# Set the filename and path for the output plot
output_filename <- "raw_and_smoothed_data_day1_100_Basis.png"
# Open a PNG device to save the plot
png(filename = output_filename, width = 800, height = 800, res = 140)
plot_raw_and_smoothed_data(data_wide$MovementAcceleration_g_d5, smoothed_data_orig_3, day = 5)
dev.off()




# Step 3: Normalize data
normalized_data <- normalize_data(data_wide)
smoothed_data <- smooth_data(normalized_data, nbasis = 25, basis_type = "fourier")
plot_raw_and_smoothed_data(normalized_data$MovementAcceleration_g_d1, smoothed_data, day = 1)

# Step 4: Create mean plots
mean_function_raw <- calculate_mean_function(normalized_data)
mean_function_df <- data.frame(t_values = normalized_data$t_values, mean_function = mean_function_raw)
mean_function_raw_1 <- calculate_mean_function(data_wide)
mean_function_df_1 <- data.frame(t_values = data_wide$t_values, mean_function = mean_function_raw_1)

# Smooth all days and create mean smoothed plot
smoothed_mean_fd <- smooth_data(mean_function_raw, nbasis = 20, basis_type = "fourier")
#plot_all_functions_with_mean_and_smoothed_mean(normalized_data, mean_function_raw, smoothed_mean_fd)
smoothed_mean_fd_1 <- smooth_data(mean_function_raw_1, nbasis = 20, basis_type = "fourier")
plot_all_functions_with_mean_and_smoothed_mean(data_wide, mean_function_raw_1, smoothed_mean_fd_1)


# Take the average over raw data of all days and smooth this function
#mean_fd <- mean.fd(smoothed_data$fd)
#plot_all_functions_with_mean(smoothed_data)
mean_fd_1 <- mean.fd(smoothed_data_orig$fd)
plot_all_functions_with_mean(smoothed_data_orig)





### split into work = non_work
work_status <- c(1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0) ##example 
setwd("/Users/nicol/Desktop/SP5/data/ID Results für Nicolas Volz")
# Step 1: Read data
data_wide <- read_and_transform_data(file_path, xlsx = FALSE)

# Step 2: Split data into work and non-work days
split_data <- split_data_by_work_status(data_wide, work_status)

work_data <- split_data$work_data
non_work_data <- split_data$non_work_data

# Continue with processing the work data
smoothed_data_work <- smooth_data(work_data, nbasis = 25, basis_type = "fourier")
plot_raw_and_smoothed_data(work_data$MovementAcceleration_g_d3, smoothed_data_work, day = 1)

# Continue with processing the non-work data
smoothed_data_non_work <- smooth_data(non_work_data, nbasis = 25, basis_type = "fourier")
plot_raw_and_smoothed_data(non_work_data$MovementAcceleration_g_d4, smoothed_data_non_work, day = 1)

# Create mean plots for work and non-work data
mean_function_work <- calculate_mean_function(work_data)
mean_function_df_work <- data.frame(t_values = work_data$t_values, mean_function = mean_function_work)

mean_function_non_work <- calculate_mean_function(non_work_data)
mean_function_df_non_work <- data.frame(t_values = non_work_data$t_values, mean_function = mean_function_non_work)

# Smooth mean functions
smoothed_mean_fd_work <- smooth_data(mean_function_work, nbasis = 25, basis_type = "fourier")
smoothed_mean_fd_non_work <- smooth_data(mean_function_non_work, nbasis = 25, basis_type = "fourier")

# Plot all functions with mean and smoothed mean for work data
plot_all_functions_with_mean_and_smoothed_mean(work_data, mean_function_work, smoothed_mean_fd_work, title = "Work Days")

# Plot all functions with mean and smoothed mean for non-work data
plot_all_functions_with_mean_and_smoothed_mean(non_work_data, mean_function_non_work, smoothed_mean_fd_non_work, title = "Non-Work Days")


plot_all_functions_with_mean(smoothed_data_work)
plot_all_functions_with_mean(smoothed_data_non_work)

plot_combined_functions_with_mean(smoothed_data_work, smoothed_data_non_work)







#################### Analysis with complete data set ###########################

setwd("/Users/nicol/Desktop/SP5/data/ID Results für Nicolas Volz/NicolasVolz_Results_IDs/")

# List all .xlsx files in the directory
file_list <- list.files(pattern = "*.xlsx")

# Read all files into a list of data frames
data_list <- lapply(file_list, function(file) read_data_file(file))


# Get the number of rows for each file
row_counts <- count_rows_in_files(data_list)

# Display the number of rows for each file
row_counts_df <- data.frame(File = file_list, RowCount = row_counts)
print(row_counts_df)

# Identify files with less than 18720 rows
files_to_keep <- row_counts_df %>%
  filter(RowCount >= 18720) %>%
  pull(File)



# Update the file_list to only include files with 18720 rows or more
file_list <- files_to_keep



### split in pain / no_pain 

extracted_ids <- sapply(file_list, extract_id)


# Standardize the extracted IDs
standardized_extracted_ids <- sapply(extracted_ids, standardize_id)
length(standardized_extracted_ids)


# Split the file list based on the extracted IDs
pain_files <- file_list[extracted_ids %in% pain_ids]
no_pain_files <- file_list[extracted_ids %in% no_pain_ids]
men_files <- file_list[extracted_ids %in% men_ids]
women_files <- file_list[extracted_ids %in% women_ids]



length(pain_files)
length(no_pain_files)


## bring in right format for the names of excel sheets


not_in_pain_or_no_pain <- file_list[!(standardized_extracted_ids %in% c(pain_ids, no_pain_ids))]




setwd("/Users/nicol/Desktop/SP5/data/ID Results für Nicolas Volz/NicolasVolz_Results_IDs/")
##################### split work / no_work 

# Run the processing for all files
work_status_list <- final_list

# Process all valid files
results_1 <- process_all_files_1(file_list, work_status_list)
results_1_pain <- process_all_files_1(pain_files, work_status_list)
results_1_no_pain <- process_all_files_1(no_pain_files, work_status_list)


setwd("/Users/nicol/Desktop/SP5/Results/")


results_1 <- Filter(Negate(is.null), results_1)
results_1_pain <- Filter(Negate(is.null), results_1_pain)
results_1_no_pain <- Filter(Negate(is.null), results_1_no_pain)


#plot_mean_and_smoothed_mean(filtered_results$`2001_u300_Results.xlsx`, "combined")
plot_mean_and_smoothed_mean(results_1$u300$work, "work")
plot_mean_and_smoothed_mean(results_1$u300$non_work, "non_work")


# Split results
split_results_1 <- split_results(results_1)
split_results_1_pain <- split_results(results_1_pain)
split_results_1_no_pain <- split_results(results_1_no_pain)


## plot all function unfiltered: 
par(mfrow = c(1,1))
output_file <- "unfiltered_data_plot.png"
# Open a PNG device to start capturing the plot
png(filename = output_file, width = 1500, height = 1500, res = 180)
# Create the plot with the specified function
plot_all_raw_and_overall_mean_functions(split_results_1$work, title = "Unfiltered work data")
# Close the PNG device to save the file
dev.off()                                  


# Filter the results for work and non-work data
filtered_results_work <- filter_results(split_results_1$work)
filtered_results_non_work <- filter_results(split_results_1$non_work)
filtered_results_pain_work <- filter_results(split_results_1_pain$work)
filtered_results_pain_non_work <- filter_results(split_results_1_pain$non_work)
filtered_results_no_pain_work <- filter_results(split_results_1_no_pain$work)
filtered_results_no_pain_non_work <- filter_results(split_results_1_no_pain$non_work)

# Create the overall mean function for work and non-work data
overall_mean_result_work <- create_overall_mean_function(filtered_results_work)
overall_mean_result_non_work <- create_overall_mean_function(filtered_results_non_work)
overall_mean_result_pain_work <- create_overall_mean_function(filtered_results_pain_work)
overall_mean_result_pain_non_work <- create_overall_mean_function(filtered_results_pain_non_work)
overall_mean_result_no_pain_work <- create_overall_mean_function(filtered_results_no_pain_work)
overall_mean_result_no_pain_non_work <- create_overall_mean_function(filtered_results_no_pain_non_work)


#### which ids when problematic
#### my mistake
# Use the modified function for each data subset
filtered_work_data <- filter_results_with_violations(split_results_1$work)
filtered_results_work <- filtered_work_data$filtered_results
excluded_work_with_violations <- filtered_work_data$exclusions

# Use the function to print each exclusion list
#print_exclusions(excluded_work_with_violations, "work")


## plot all function unfiltered: 
par(mfrow = c(1,1))
output_file <- "Filtered_data_plot.png"
# Open a PNG device to start capturing the plot
png(filename = output_file, width = 1500, height = 1500, res = 180)
# Create the plot with the specified function
plot_all_raw_and_overall_mean_functions(filtered_results_work, title = "Filtered work data")
# Close the PNG device to save the file
dev.off()                                  














#############-----build FoSR Model-----#########################################


# work day model 
### design Matrix X 

## info on covariates

# Define evaluation grid
evaluation_grid <- seq(0, 1440, by = 1)  # 1441 points



setwd("/Users/nicol/Desktop/SP5/Results/")


to_merge = data.frame(id = merged_df_3$id, age=merged_df_3$age, gender=merged_df_3$gender)
to_merge = na.omit(to_merge)
### work day model 
work_day_model = do_regression(
  pain_data = filtered_results_pain_work,
  nopain_data = filtered_results_no_pain_work,
  n_pain = length(filtered_results_pain_work),
  n_nopain = length(filtered_results_no_pain_work),
  id = c(names(filtered_results_pain_work), names(filtered_results_no_pain_work)),
  to_merge = to_merge,
  evaluation_grid = evaluation_grid,
  model_name = "Work_Days", 
  plot = TRUE
)
work_day_model[2]



### non work day model 
non_work_day_model = do_regression(
  pain_data = filtered_results_pain_non_work,
  nopain_data = filtered_results_no_pain_non_work,
  n_pain = length(filtered_results_pain_non_work),
  n_nopain = length(filtered_results_no_pain_non_work),
  id = c(names(filtered_results_pain_non_work), names(filtered_results_no_pain_non_work)),
  to_merge = to_merge,
  evaluation_grid = evaluation_grid,
  model_name = "Non_Work_Days", 
  plot = TRUE
)
non_work_day_model[2]



#### pain characterized in different ways below################################

'
###### kroff 
to_merge_kroff = data.frame(id = merged_df_3$id, age=merged_df_3$age, gender=merged_df_3$gender, kroff = merged_df_3$intensitaet)
### work day model_kroff
work_day_model_kroff = do_regression(
  pain_data = filtered_results_pain_work,
  nopain_data = filtered_results_no_pain_work,
  n_pain = length(filtered_results_pain_work),
  n_nopain = length(filtered_results_no_pain_work),
  id = c(names(filtered_results_pain_work), names(filtered_results_no_pain_work)),
  to_merge = to_merge_kroff,
  evaluation_grid = evaluation_grid,
  model_name = "Work_Days_kroff",
  kroff = TRUE
)
work_day_model_kroff[2]


### non work day model kroff
non_work_day_model_kroff = do_regression(
  pain_data = filtered_results_pain_non_work,
  nopain_data = filtered_results_no_pain_non_work,
  n_pain = length(filtered_results_pain_non_work),
  n_nopain = length(filtered_results_no_pain_non_work),
  id = c(names(filtered_results_pain_non_work), names(filtered_results_no_pain_non_work)),
  to_merge = to_merge_kroff,
  evaluation_grid = evaluation_grid,
  model_name = "Non_Work_Days_kroff", 
  kroff = TRUE
)
non_work_day_model[2]



###### beeintraechtigung 
to_merge_kroff = data.frame(id = merged_df_3$id, age=merged_df_3$age, gender=merged_df_3$gender, kroff = merged_df_3$beeintraechtigung)
### work day model_beeintraechtigung
work_day_model_kroff = do_regression(
  pain_data = filtered_results_pain_work,
  nopain_data = filtered_results_no_pain_work,
  n_pain = length(filtered_results_pain_work),
  n_nopain = length(filtered_results_no_pain_work),
  id = c(names(filtered_results_pain_work), names(filtered_results_no_pain_work)),
  to_merge = to_merge_kroff,
  evaluation_grid = evaluation_grid,
  model_name = "Work_Days_beein",
  kroff = TRUE
)
work_day_model_kroff[2]


### non work day model beeintraechtigung
non_work_day_model = do_regression(
  pain_data = filtered_results_pain_non_work,
  nopain_data = filtered_results_no_pain_non_work,
  n_pain = length(filtered_results_pain_non_work),
  n_nopain = length(filtered_results_no_pain_non_work),
  id = c(names(filtered_results_pain_non_work), names(filtered_results_no_pain_non_work)),
  to_merge = to_merge_kroff,
  evaluation_grid = evaluation_grid,
  model_name = "Non_Work_Days_beein", 
  kroff = TRUE
)
non_work_day_model[2]
'




################################################################################
############################################### Differance Model################ 
################################################################################


fd_work = work_day_model[[1]]
fd_nonwork = non_work_day_model[[1]]
mean_fd_work = mean.fd(fd_work)
mean_fd_nonwork = mean.fd(fd_nonwork)


# Extract subjects with at least two work and two non-work days
subject_ids <- unique(c(names(filtered_results_work), names(filtered_results_non_work)))

valid_subjects <- sapply(subject_ids, function(id) {
  work_days <- sum(work_status_list[[id]] == 1)
  non_work_days <- sum(work_status_list[[id]] == 0)
  return(work_days >= 2 && non_work_days >= 2)
})

# Get the IDs of valid subjects
valid_subject_ids <- names(valid_subjects)[valid_subjects]


## split into pain / np pain

# Split the file list based on the extracted ID
valid_subject_ids_pain <- valid_subject_ids[valid_subject_ids %in% pain_ids]
valid_subject_ids_nopain <- valid_subject_ids[valid_subject_ids %in% no_pain_ids]



split_dif_ids = c(valid_subject_ids_pain, valid_subject_ids_nopain)


# Initialize a list to store the difference functions
difference_functions <- list()

# Loop through valid subjects to compute difference functions
for (id in split_dif_ids) {
  work_mean_function <- (filtered_results_work[[id]]$mean_function_raw)
  non_work_mean_function <- (filtered_results_non_work[[id]]$mean_function_raw)
  difference_function <- work_mean_function - non_work_mean_function
  difference_functions[[id]] <- difference_function
}

plot(seq(1,1440,by=1),difference_functions$u098, type = "l")


# Initialize an empty list to store the coefficients
smoothed_difference_coefs <- list()

# Loop through valid subjects to compute and smooth difference functions
for (id in split_dif_ids) {
  # Compute mean functions for work and non-work days
  work_mean_function <- filtered_results_work[[id]]$mean_function_raw
  non_work_mean_function <- filtered_results_non_work[[id]]$mean_function_raw
  
  # Compute the difference function
  difference_function <- work_mean_function - non_work_mean_function
  
  # Smooth the difference function using the smooth_data function
  smoothed_difference_result <- smooth_data(difference_function, nbasis = 25, basis_type = "fourier")
  
  # Store the coefficients of the smoothed function
  smoothed_difference_coefs[[id]] <- smoothed_difference_result$fd$coefs
}

# Combine the coefficients into one big matrix
# Each column represents the smoothed coefficients for a subject
coef_matrix_dif <- do.call(cbind, smoothed_difference_coefs)


# Check the dimensions of the combined matrix to ensure it's correct
print(dim(coef_matrix_dif))


# Now create a B-spline basis object for the fd object
rangeval <- c(0, 1440)  # The range of your data
nbasis <- 25  # Number of B-spline basis functions

# Create a B-spline basis object
basisobj <- create.fourier.basis(rangeval, nbasis)

# Create the fd object using the combined coefficients and the basis object
fd_obj_dif <- fd(coef = coef_matrix_dif, basisobj = basisobj)



results_intercept_only <- do_intercept_only_regression(coef_matrix_dif, evaluation_grid, "Difference_Model")






################### OVERAL MEASURE BY STEPS#####################################

setwd("/Users/nicol/Desktop/SP5/data/ID Results für Nicolas Volz/NicolasVolz_Results_IDs/")

############################ step counts 
results_steps <- process_all_files_steps(file_list, work_status_list)
results_pain_steps <- process_all_files_steps(pain_files, work_status_list)
results_no_pain_steps <- process_all_files_steps(no_pain_files, work_status_list)
results_men <- process_all_files_steps(men_files, work_status_list)
results_women <- process_all_files_steps(women_files, work_status_list)

results_means <- calculate_mean_work_nonwork_steps(results_steps, work_status_list)
results_means_pain <- calculate_mean_work_nonwork_steps(results_pain_steps, work_status_list)
results_means_no_pain <- calculate_mean_work_nonwork_steps(results_no_pain_steps, work_status_list)
results_means_men <- calculate_mean_work_nonwork_steps(results_men, work_status_list)
results_means_women <- calculate_mean_work_nonwork_steps(results_women, work_status_list)


results_steps$mean_steps <- rowMeans(results_steps[,-1], na.rm = TRUE)  # Exclude the first column (ID)
results_pain_steps$mean_steps <- rowMeans(results_pain_steps[,-1], na.rm = TRUE)  # Exclude the first column (ID)
results_no_pain_steps$mean_steps <- rowMeans(results_no_pain_steps[,-1], na.rm = TRUE)  # Exclude the first column (ID)
results_men$mean_steps <- rowMeans(results_men[,-1], na.rm = TRUE)  # Exclude the first column (ID)
results_women$mean_steps <- rowMeans(results_women[,-1], na.rm = TRUE)  # Exclude the first column (ID)



setwd("/Users/nicol/Desktop/SP5/Results/")

# Initialize an empty data frame to hold all the results
all_stats <- data.frame()

# Example 1: results_means$work_mean vs. results_means$non_work_mean
result_work_vs_non_work <- plot_boxplots_comparison(results_means$work_mean, results_means$non_work_mean, 
                                                    group1_label = "Work Days", group2_label = "Non-Work Days", 
                                                    title = "")

# Save plot and get stats
all_stats <- rbind(all_stats, save_plot_and_return_stats(result_work_vs_non_work, "plot_work_vs_non_work.png"))

# Example 2: results_means_pain$work_mean vs. results_means_no_pain$work_mean
result_pain_vs_no_pain <- plot_boxplots_comparison(results_means_pain$work_mean, results_means_no_pain$work_mean, 
                                                   group1_label = "Pain", group2_label = "No Pain", 
                                                   title = "")

# Save plot and get stats
all_stats <- rbind(all_stats, save_plot_and_return_stats(result_pain_vs_no_pain, "plot_pain_vs_no_pain_work.png"))

# Example 3: results_means_pain$non_work_mean vs. results_means_no_pain$non_work_mean
result_pain_vs_no_pain_non_work <- plot_boxplots_comparison(
  results_means_pain$non_work_mean, results_means_no_pain$non_work_mean, 
  group1_label = "Pain", group2_label = "No Pain", 
  title = ""
)

# Save plot and get stats
all_stats <- rbind(all_stats, save_plot_and_return_stats(result_pain_vs_no_pain_non_work, "plot_pain_vs_no_pain_non_work.png"))

# Example 4: results_means_men$non_work_mean vs. results_means_women$non_work_mean
result_men_vs_women_non_work <- plot_boxplots_comparison(
  results_means_men$non_work_mean, results_means_women$non_work_mean, 
  group1_label = "Men", group2_label = "Women", 
  title = ""
)

# Save plot and get stats
all_stats <- rbind(all_stats, save_plot_and_return_stats(result_men_vs_women_non_work, "plot_men_vs_women_non_work.png"))

# Example 5: results_means_men$work_mean vs. results_means_women$work_mean
result_men_vs_women_work <- plot_boxplots_comparison(
  results_means_men$work_mean, results_means_women$work_mean, 
  group1_label = "Men", group2_label = "Women", 
  title = ""
)

# Save plot and get stats
all_stats <- rbind(all_stats, save_plot_and_return_stats(result_men_vs_women_work, "plot_men_vs_women_work.png"))

# Export the results to a CSV file
write.csv(all_stats, "comparison_results.csv", row.names = FALSE)


