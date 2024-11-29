
#### packages 
library(dplyr)



## fix problem that some ids are e.g. 080 and some are 80
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


setwd("data info")

################################################################################
##############################---read raw data-----#############################
################################################################################

### chronic pain? 
T0_raw <- read.csv(file = "Daten_SP5_T0_aufbereitet_240517.csv", header = TRUE, sep = "," , dec = ",")
### fix leading zero Problem: id: 83 to id:083
T0_raw$sp5_id <- sapply(T0_raw$sp5_id, format_id)
T0 <- data.frame(id = T0_raw$sp5_id, pain = T0_raw$sp5_chrpain)



dayly_pain_raw <- read.csv(file = "Daily_data_Nicolas_240701.csv", header = TRUE, sep = ";" , dec = ",")
dayly_pain <- distinct(data.frame(id = dayly_pain_raw$sp5_id, id2 = dayly_pain_raw$ID))


### personal info
Covariates_raw <- read.csv(file = "T0_Daten_Nicolas_240628.csv", header = TRUE, sep = ";" , dec = ",")
Covariates_raw$sp5_id <- sapply(Covariates_raw$sp5_id, format_id)
Covariates <- data.frame(id = Covariates_raw$sp5_id, age = Covariates_raw$sp5_age, 
                         gender = Covariates_raw$sp5_gender, intensitaet = Covariates_raw$sp5_korff_schmerzintensitaet,
                         beeintraechtigung = Covariates_raw$sp5_korff_beeintraechtigung, dis_score = Covariates_raw$sp5_korff_dis_score)
Covariates$id = sapply(Covariates$id, format_id)


### pain eps per day
Pain_eps_raw <- read.csv(file = "Data_number_of_pain_episodes_per_ID_240701.csv", header = TRUE, sep = ";" , dec = ",")
Pain_eps <- data.frame(id2 = Pain_eps_raw$ID, pain_eps = Pain_eps_raw$pain_epis_number)



### work or no work day. problematic because very few work days
### the data is split here mo-fr -> work. sa-so -> free 
Work_raw <- read.csv(file = "Daily_data_mp5_only_Nicolas_240701.csv", header = TRUE, sep = ";" , dec = ",", na.strings = "NA")
Work_raw$day <- as.Date(Work_raw$date, format = "%d.%m.%Y")
Work_raw$day_of_week <- weekdays(Work_raw$day)
Work_raw$workday <- ifelse(weekdays(Work_raw$day) %in% c("Samstag", "Sonntag"), 0, 1)

#Work <- data.frame(id2 = Work_raw$ID, date = Work_raw$date, arbeit = Work_raw$arbeit1)
Work <- data.frame(id2 = Work_raw$ID, date = Work_raw$date, arbeit = Work_raw$workday)
### ACHTUNG: z.B. a_11 -> a_011 verändert in T0_Daten_Nicolas_240628.csv (no real ids) 




####################-- combine data---##########################################
#### 
merged_df_1 <- merge(dayly_pain, T0, by = "id", all = TRUE)
merged_df_2 <- merge(merged_df_1, Covariates, by = "id", all = TRUE)
merged_df_3 <- merge(merged_df_2, Pain_eps, by = "id2", all = TRUE)
merged_df_4 <- merge(merged_df_3, Work, by = "id2", all = TRUE)
merged_df_4$arbeit[is.na(merged_df_4$arbeit)] <- 0 






############-------- define pain groups-------################################### 
## identify pain group
mean(merged_df_4$intensitaet[merged_df_4$pain==1],na.rm = TRUE)
mean(merged_df_4$intensitaet[merged_df_4$pain==2],na.rm = TRUE)

mean(merged_df_4$beeintraechtigung[merged_df_4$pain==1],na.rm = TRUE)
mean(merged_df_4$beeintraechtigung[merged_df_4$pain==2],na.rm = TRUE)

mean(merged_df_4$dis_score[merged_df_4$pain==1],na.rm = TRUE)
mean(merged_df_4$dis_score[merged_df_4$pain==2],na.rm = TRUE)

mean(merged_df_4$pain_eps[merged_df_4$pain==1],na.rm = TRUE)
mean(merged_df_4$pain_eps[merged_df_4$pain==2],na.rm = TRUE)

## 1 -> pain 
## 2 -> no pain 



## pain ids anhand chronizität
## pain and no pain grup
merged_pain_chron <- merged_df_3[merged_df_3$pain == 1 & !is.na(merged_df_3$pain), ]
merged_no_pain_chron <- merged_df_3[merged_df_3$pain == 2 & !is.na(merged_df_3$pain), ]


pain_ids_chron = merged_pain_chron$id
no_pain_ids_chron = merged_no_pain_chron$id

## pain ids anhand von intensität (kroff_40)
## pain and no pain grup
merged_pain_kroff_40 <- merged_df_3[merged_df_3$intensitaet >=40 & !is.na(merged_df_3$intensitaet), ]
merged_no_pain_kroff_40 <- merged_df_3[merged_df_3$intensitaet < 40 & !is.na(merged_df_3$intensitaet), ]


pain_ids_kroff_40 = merged_pain_kroff_40$id
no_pain_ids_kroff_40 = merged_no_pain_kroff_40$id

## pain ids anhand von intensität (kroff_50)
## pain and no pain grup
merged_pain_kroff_50 <- merged_df_3[merged_df_3$intensitaet >=50 & !is.na(merged_df_3$intensitaet), ]
merged_no_pain_kroff_50 <- merged_df_3[merged_df_3$intensitaet < 50 & !is.na(merged_df_3$intensitaet), ]


pain_ids_kroff_50 = merged_pain_kroff_50$id
no_pain_ids_kroff_50 = merged_no_pain_kroff_50$id




## pain ids anhand von beeinträchtigng 
## pain and no pain grup
merged_pain_beeintraechtigung  <- merged_df_3[merged_df_3$beeintraechtigung  >=30 & !is.na(merged_df_3$beeintraechtigung ), ]
merged_no_pain_beeintraechtigung  <- merged_df_3[merged_df_3$beeintraechtigung  < 30 & !is.na(merged_df_3$beeintraechtigung ), ]


pain_ids_beeintraechtigung  = merged_pain_beeintraechtigung $id
no_pain_ids_beeintraechtigung  = merged_no_pain_beeintraechtigung $id


## pain ids anhand von auffaellig 
## pain and no pain grup
merged_pain_auffaellig   <- merged_df_3[merged_df_3$beeintraechtigung  >=30 & !is.na(merged_df_3$beeintraechtigung ), ]
merged_no_pain_auffaellig   <- merged_df_3[merged_df_3$beeintraechtigung  < 30 & !is.na(merged_df_3$beeintraechtigung ), ]


pain_ids_auffaellig   = merged_pain_beeintraechtigung $id
no_pain_ids_auffaellig   = merged_no_pain_beeintraechtigung $id




## pain ids anhand von beeinträchtigng 
## pain and no pain grup
merged_pain_extreme  <- merged_df_3[merged_df_3$intensitaet  >=50 & !is.na(merged_df_3$intensitaet ), ]
merged_no_pain_extreme  <- merged_df_3[merged_df_3$intensitaet  < 5 & !is.na(merged_df_3$intensitaet), ]


pain_ids_extreme  = merged_pain_extreme $id
no_pain_ids_extreme  = merged_no_pain_extreme $id



###############---- men and women------#########################################


men = merged_df_3[merged_df_3$gender==1& !is.na(merged_df_3$gender ),]
women = merged_df_3[merged_df_3$gender==2 &!is.na(merged_df_3$gender ),]

men_ids = men$id
women_ids = women$id



################# split into work no_work#######################################


# Group by id and summarize arbeit into a list with names
work_status_list <- merged_df_4 %>%
  group_by(id) %>%
  summarize(work_status = list(arbeit)) %>%
  deframe()



## only keep for ones that are longer than 13 days. only look at first 13 days
## because day 14 often incomplete 
filtered_list <- work_status_list[sapply(work_status_list, length) >= 13]
final_list <- lapply(filtered_list, function(x) x[-14])


######################### checks ###############################################

### correctly meged? 
id_1 = (unique(merged_df_4$id[merged_df_4$pain==2]))
id_2 = (unique(merged_df_4$id2[merged_df_4$pain==2]))

id_11 = (unique(merged_df_1$id[merged_df_1$pain==2]))
id_22 = (unique(merged_df_1$id2[merged_df_1$pain==2]))

setequal(id_1, id_11) # true

id_2_problem = merged_df_1$id2[merged_df_1$id  %in%  c(id_1)]





### mean and median for work days and non work days 
# Function to calculate workdays (1's) and non-workdays (0's) for each person
workdays_count <- sapply(work_status_list, function(x) sum(x == 1))
non_workdays_count <- sapply(work_status_list, function(x) sum(x == 0))

# Calculate the mean and median number of workdays
mean_workdays <- mean(workdays_count)
median_workdays <- median(workdays_count)

# Calculate the mean and median number of non-workdays
mean_non_workdays <- mean(non_workdays_count)
median_non_workdays <- median(non_workdays_count)

# Print the results
#cat("Mean workdays:", mean_workdays, "\n")
#cat("Median workdays:", median_workdays, "\n")
#cat("Mean non-workdays:", mean_non_workdays, "\n")
#cat("Median non-workdays:", median_non_workdays, "\n")




