# Load required library
library(dplyr)
library(tidyr)

source('basic_function.R')


##### ECC ####

energy_function <- function(eval='smooth',python_df=TRUE){
  setwd("./energy")

  
  
  #########################################################
  #########################################################
  ## full categories
  #########################################################
  #########################################################
  df <- read.csv("world_energy.csv")
  
  energy_cols <- c("coalcons_ej", "gascons_ej", "oilcons_ej", 
                   "nuclear_ej", "hydro_ej", "renewables_ej")
  
  # 3. Calculate percentages
  df_pct <- df %>%
    filter(!grepl("^Total|^Rest of|^Other", Country)) %>%
    
    select(Country, Year, all_of(energy_cols)) %>%
    mutate(total_energy = rowSums(select(., all_of(energy_cols)), na.rm = TRUE)) %>%
    mutate(across(all_of(energy_cols), ~ (. / total_energy))) %>%
    select(-total_energy)
  
  
  country_list <- split(df_pct, df_pct$Country)
  
  country_matrices <- lapply(country_list, function(country_df) {
    mat <- as.matrix(country_df[, energy_cols])
    rownames(mat) <- country_df$Year
    
    return(mat)
  })
  
  country_matrices_clean <- Filter(function(mat) !any(is.na(mat)), country_matrices)
  
  years_list <- lapply(country_matrices_clean, rownames)
  
  # Use Reduce and intersect to find the overlapping years
  common_years <- Reduce(intersect, years_list)
  
  
  # --- Outputs ---
  
  cat("Original number of countries:", length(country_matrices), "\n")
  cat("Clean number of countries (no NAs):", length(country_matrices_clean), "\n\n")
  
  cat("Common years across clean countries:\n")
  # Convert the rownames back to numeric for a clean sorted print
  print(sort(as.numeric(common_years)))
  
  
  #####################################
  ## dissimilarity matrix
  full_list <- country_matrices_clean
  dissimilarity_matrix <- matrix(0,length(full_list),length(full_list))
  ctr_name <- names(country_matrices_clean)
  row.names(dissimilarity_matrix) = colnames(dissimilarity_matrix) = ctr_name
  
  naive_dissimilarity_matrix <- dissimilarity_matrix
  SPD_dissimilarity_matrix <- dissimilarity_matrix  ##Spherical Procrustes distance
  ESD_dissimilarity_matrix <- dissimilarity_matrix ##Euclidean shape distance
  ED_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean distance
  EPD_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean Procrustes distance
  ##################################################
  
  if(eval=='smooth'){
    ## smooth
    smooth_list <- list()
    
    for (i in 1:(length(full_list))) {
      sub_matrix_i <- sqrt(full_list[[i]])
      
      sample_size <- nrow(sub_matrix_i)
      x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1)
      
      trend_component <- LocSpheReg(xin = as.vector(x_input), yin = sub_matrix_i)
      bw_current <- trend_component$optns$bw 
      print(bw_current)
      sub_matrix_smooth <- trend_component$yout
      row.names(sub_matrix_smooth) <- row.names(sub_matrix_i)
      smooth_list[[i]] <- sub_matrix_smooth
      print(i)
    }
    names(smooth_list) <- ctr_name
    
  } else{
    
    ## no smooth
    smooth_list <- list()
    
    for (i in 1:(length(full_list))) {
      sub_matrix_i <- sqrt(full_list[[i]])
      smooth_list[[i]] <- sub_matrix_i
      print(i)
    }
    names(smooth_list) <- ctr_name
    
  }
  
  
  ##################################################
  
  mean_matrix <- matrix(NA,nrow=length(full_list),ncol=ncol(smooth_list[[i]]))
  for (i in 1:(length(full_list))) {
    sub_matrix_i <- smooth_list[[i]]
    mean_matrix[i,] <- intrinsic_mean_sphere(sub_matrix_i)
    print(i)
  }
  common_mean <- intrinsic_mean_sphere(mean_matrix)
  
  for (i in 1:(length(full_list)-1)) {
    for (j in (i+1):length(full_list)) {
      sub_matrix_i <- smooth_list[[i]]
      sub_matrix_j <- smooth_list[[j]]
      
      ##################################################
      ##shape
      sub_matrix_i_transport <- transport_shape_function(sub_matrix_i,common_mean)[[1]]
      sub_matrix_j_transport <- transport_shape_function(sub_matrix_j,common_mean)[[1]]
      dis_measure <- shape_dissimilarity_function_1D(sub_matrix_i_transport,sub_matrix_j_transport)
      dissimilarity_matrix[i,j] <- dis_measure
      dissimilarity_matrix[j,i] <- dis_measure
      
      ##naive great circle distance
      dis_measure <- sphere_ts_dissimilarity_function(sub_matrix_i,sub_matrix_j)
      naive_dissimilarity_matrix[i,j] <- dis_measure
      naive_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Spherical Procrustes distance
      dis_measure <- sphere_procrustes_distance(sub_matrix_i,sub_matrix_j)
      SPD_dissimilarity_matrix[i,j] <- dis_measure
      SPD_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Euclidean shape distance
      dis_measure <- multi_function_dissimilarity_function(sub_matrix_i,sub_matrix_j)
      ESD_dissimilarity_matrix[i,j] <- dis_measure
      ESD_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Euclidean distance
      dis_measure <- function_shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
      ED_dissimilarity_matrix[i,j] <- dis_measure
      ED_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Euclidean Procrustes distance
      dis_measure <- procrustes_distance(sub_matrix_i,sub_matrix_j)
      EPD_dissimilarity_matrix[i,j] <- dis_measure
      EPD_dissimilarity_matrix[j,i] <- dis_measure
    }
    print(i)
  }
  
  ##### classification
  country_labels <- data.frame(
    ctr_name = ctr_name,
    oecd_status = case_match(
      ctr_name,
      c("Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia", "Czech Republic", "Denmark", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Israel", "Italy", "Japan", "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway", "Poland", "Portugal", "Slovakia", "South Korea", "Spain", "Sweden", "Switzerland", "Turkiye", "United Kingdom", "US") ~ "OECD",
      .default = "Non-OECD"
    )
  )
  
  
  target_label <- "oecd_status" 
  # ---------------------------------------------------------
  # 1. Extract All Items and Map Labels
  # ---------------------------------------------------------
  all_items <- rownames(dissimilarity_matrix)

  all_labels <- country_labels[[target_label]][match(all_items, country_labels$ctr_name)]
  names(all_labels) <- all_items
  
  # ---------------------------------------------------------
  # 2. Create a Train / Test Split
  # ---------------------------------------------------------
  set.seed(123) 
  
  n_total <- length(all_items)
  train_size <- floor(0.75 * n_total)
  
  train_items <- sample(all_items, size = train_size)
  test_items <- setdiff(all_items, train_items)
  
  train_labels <- all_labels[train_items]
  true_test_labels <- all_labels[test_items]  
  
  # ---------------------------------------------------------
  # 3. kNN  for Distance Matrices
  # ---------------------------------------------------------
  result_matrix <- matrix(0,nrow = 2,ncol = 6)
  ##SSD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,1] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(dissimilarity_matrix,
                                                all_items,
                                                all_labels,
                                                train_items,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,1] <- round(results[[4]] * 100, 2)


  ##SD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = naive_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,2] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(naive_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,2] <- round(results[[4]] * 100, 2)
  
  ##SPD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = SPD_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,3] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(SPD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,3] <- round(results[[4]] * 100, 2)
  
  ##ESD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = ESD_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,4] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(ESD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 

  result_matrix[2,4] <- round(results[[4]] * 100, 2)
  
  
  ##ED
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = ED_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,5] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(ED_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,5] <- round(results[[4]] * 100, 2)
  
  ##EPD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = EPD_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,6] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(EPD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,6] <- round(results[[4]] * 100, 2)
  
  
  ## prepare for Python 
  if (python_df==TRUE){
    names(smooth_list) <- country_labels$oecd_status[match(names(smooth_list), country_labels$ctr_name)]
    save(smooth_list,file = 'matrix_input.RData')

    dir.create("matrix_export", showWarnings = FALSE)

    for (i in seq_along(smooth_list)) {
      label_name <- names(smooth_list)[i]
      file_name <- paste0("matrix_export/", label_name, "_", i, ".csv")

      write.table(smooth_list[[i]], file = file_name, sep = ",", 
                  row.names = FALSE, col.names = FALSE)
    } 
  }
  return(result_matrix)
}

##Table 1 part 1 for ECC
energy_function(eval = 'smooth')

##Table 4 part 1 for ECC
energy_function(eval = 'non-smooth')

##### LTDC #####
# First downlowd the LDTC for male via https://www.mortality.org and put it in the folder 'HMD'
HMD_function <- function(eval='smooth',python_df=TRUE){
  setwd("./HMD")
  
  excluded_ctr <- c('BEL', 'UKR', 'RUS', 'SVN', 'NZL_NM', 'NZL_MA', 'KOR', 'ISR', 'HRV', 'GRC', 'CHL')
  # male
  ages = 0:110
  male_table <- read.table("mltper_1x1.txt", skip=2, header= TRUE)
  male_multicontry_list <- list()
  year_multicontry_list <- list()
  
  for (ctr in unique(male_table$PopName)) {
    if (!(ctr %in% excluded_ctr)) {
      current_table <- male_table[male_table$PopName==ctr,]
      years <- unique(current_table$Year)
      year_multicontry_list[[ctr]] <- years
      ctr_table <- t(matrix(current_table[,"qx"], 111, ))
      male_pop <- matrix(NA, nrow(ctr_table), ncol(ctr_table))
      for(ij in 1:nrow(ctr_table))
      {
        # set radix
        start_pop = 10^5
        for(ik in 1:ncol(ctr_table))
        {
          raw_val <- unlist(ctr_table[ij, ik])
          num_val <- as.numeric(raw_val)
          
          # Report if it turned into NA but wasn't originally empty
          if (is.na(num_val) && !is.na(raw_val) && raw_val != "") {
            message(sprintf("Coercion Issue: Row [%s], Col [%s] | Value: '%s' | Country: %s", 
                            ij, ik, raw_val, ctr))
          }
          male_pop[ij,ik] = as.numeric(unlist(ctr_table[ij,ik])) * start_pop
          start_pop = start_pop - male_pop[ij,ik]
          
        }
      }
      rownames(male_pop) = years
      colnames(male_pop) = ages
      
      all(round(rowSums(male_pop),9) == 10^5) # TRUE
      male_multicontry_list[[ctr]] <- male_pop / 10^5
    }
  }
  
  common_years <- Reduce(intersect, year_multicontry_list)
  
  
  #####################################
  ## dissimilarity matrix
  full_list <- male_multicontry_list
  dissimilarity_matrix <- matrix(0,length(full_list),length(full_list))
  ctr_name <- setdiff(unique(male_table$PopName), excluded_ctr)
  row.names(dissimilarity_matrix) = colnames(dissimilarity_matrix) = ctr_name
  
  naive_dissimilarity_matrix <- dissimilarity_matrix
  SPD_dissimilarity_matrix <- dissimilarity_matrix  ##Spherical Procrustes distance
  ESD_dissimilarity_matrix <- dissimilarity_matrix ##Euclidean shape distance
  ED_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean distance
  EPD_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean Procrustes distance
  ##################################################
  
  if(eval=='smooth'){
    smooth_list <- list()
    
    for (i in 1:(length(full_list))) {
      years_char <- as.character(common_years)
      sub_matrix_i <- sqrt(full_list[[i]][years_char, ])
      
      sample_size <- nrow(sub_matrix_i)
      x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1)
      
      trend_component <- LocSpheReg(xin = as.vector(x_input), yin = sub_matrix_i)
      bw_current <- trend_component$optns$bw 
      print(bw_current)
      sub_matrix_smooth <- trend_component$yout
      row.names(sub_matrix_smooth) <- years_char
      smooth_list[[i]] <- sub_matrix_smooth
      print(i)
    }
    names(smooth_list) <- ctr_name
  } else{
    ## no smooth
    smooth_list <- list()
    
    for (i in 1:(length(full_list))) {
      years_char <- as.character(common_years)
      sub_matrix_i <- sqrt(full_list[[i]][years_char, ])
      smooth_list[[i]] <- sub_matrix_i
      print(i)
    }
    names(smooth_list) <- ctr_name
    
  }
  
  
  mean_matrix <- matrix(NA,nrow=length(full_list),ncol=ncol(smooth_list[[i]]))
  for (i in 1:(length(full_list))) {
    sub_matrix_i <- smooth_list[[i]]
    mean_matrix[i,] <- intrinsic_mean_sphere(sub_matrix_i)
    print(i)
  }
  common_mean <- intrinsic_mean_sphere(mean_matrix)
  
  for (i in 1:(length(full_list)-1)) {
    for (j in (i+1):length(full_list)) {
      
      ## smooth
      sub_matrix_i <- smooth_list[[i]]
      sub_matrix_j <- smooth_list[[j]]
      sub_matrix_i_transport <- transport_shape_function(sub_matrix_i,common_mean)[[1]]
      sub_matrix_j_transport <- transport_shape_function(sub_matrix_j,common_mean)[[1]]
      dis_measure <- shape_dissimilarity_function_1D(sub_matrix_i_transport,sub_matrix_j_transport)
      
      dissimilarity_matrix[i,j] <- dis_measure
      dissimilarity_matrix[j,i] <- dis_measure
      
      ##naive great circle distance
      dis_measure <- sphere_ts_dissimilarity_function(sub_matrix_i,sub_matrix_j)
      naive_dissimilarity_matrix[i,j] <- dis_measure
      naive_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Spherical Procrustes distance
      dis_measure <- sphere_procrustes_distance(sub_matrix_i,sub_matrix_j)
      SPD_dissimilarity_matrix[i,j] <- dis_measure
      SPD_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Euclidean shape distance
      dis_measure <- multi_function_dissimilarity_function(sub_matrix_i,sub_matrix_j)
      ESD_dissimilarity_matrix[i,j] <- dis_measure
      ESD_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Euclidean distance
      dis_measure <- function_shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
      ED_dissimilarity_matrix[i,j] <- dis_measure
      ED_dissimilarity_matrix[j,i] <- dis_measure
      ##Euclidean Procrustes distance
      dis_measure <- procrustes_distance(sub_matrix_i,sub_matrix_j)
      EPD_dissimilarity_matrix[i,j] <- dis_measure
      EPD_dissimilarity_matrix[j,i] <- dis_measure
    }
    print(i)
  }
  
  ##### classification
  
  country_labels <- data.frame(
    ctr_name = ctr_name,
    historical_bloc = case_match(
      ctr_name,
      c("BGR", "BLR", "CZE", "DEUTE", "EST", "HUN", "LTU", "LVA", "POL", "SVK") ~ "Post-Communist",
      .default = "Other Developed"
    )
  )
  
  target_label <- "historical_bloc" 
  # ---------------------------------------------------------
  # 1. Extract All Items and Map Labels
  # ---------------------------------------------------------
  all_items <- rownames(dissimilarity_matrix)
  
  base_country_codes <- sub("_(male|female)$", "", all_items)

  all_labels <- country_labels[[target_label]][match(base_country_codes, country_labels$ctr_name)]
  names(all_labels) <- all_items
  
  # ---------------------------------------------------------
  # 2. Create a Train / Test Split
  # ---------------------------------------------------------
  set.seed(123) 
  developed_items <- names(all_labels)[all_labels == "Other Developed"]
  post_communist_items <- names(all_labels)[all_labels == "Post-Communist"]

  train_developed <- sample(developed_items, size = floor(0.75 * length(developed_items)))

  train_post_communist <- sample(post_communist_items, size = floor(0.75 * length(post_communist_items)))

  train_items <- c(train_developed, train_post_communist)

  test_items <- setdiff(all_items, train_items)

  train_labels <- all_labels[train_items]
  true_test_labels <- all_labels[test_items]
  
  # ---------------------------------------------------------
  # 3. kNN  for Distance Matrices
  # ---------------------------------------------------------
  result_matrix <- matrix(0,nrow = 2,ncol = 6)
  ##SSD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,1] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(dissimilarity_matrix,
                                                all_items,
                                                all_labels,
                                                train_items,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,1] <- round(results[[4]] * 100, 2)
  
  
  ##SD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = naive_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,2] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(naive_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,2] <- round(results[[4]] * 100, 2)
  
  ##SPD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = SPD_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,3] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(SPD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,3] <- round(results[[4]] * 100, 2)
  
  ##ESD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = ESD_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,4] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(ESD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,4] <- round(results[[4]] * 100, 2)
  
  
  ##ED
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = ED_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,5] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(ED_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,5] <- round(results[[4]] * 100, 2)
  
  ##EPD
  results <- repeated_stratified_evaluation(
    dissimilarity_matrix = EPD_dissimilarity_matrix,
    all_labels = all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123
  )
  result_matrix[1,6] <- round(results$mean_accuracy * 100, 2)
  
  results <- classification_evaluation_function(EPD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[2,6] <- round(results[[4]] * 100, 2)
  
  
  ## prepare for Python 
  if (python_df==TRUE){
    names(smooth_list) <- country_labels$oecd_status[match(names(smooth_list), country_labels$ctr_name)]
    save(smooth_list,file = 'matrix_input.RData')
    
    
    # Create a folder to hold the data
    dir.create("matrix_export", showWarnings = FALSE)
    
    # Loop through the list and save each matrix
    for (i in seq_along(smooth_list)) {
      # Get the name (e.g., "Non-OECD", "OECD")
      label_name <- names(smooth_list)[i]
      
      # Create a unique filename: e.g., matrix_export/Non-OECD_1.csv
      file_name <- paste0("matrix_export/", label_name, "_", i, ".csv")
      
      # Write to CSV, dropping row/column names to keep it clean
      write.table(smooth_list[[i]], file = file_name, sep = ",", 
                  row.names = FALSE, col.names = FALSE)
    } 
  }
  return(result_matrix)
}

##Table 1 part 1 for LTDC
HMD_function(eval = 'smooth')

##Table 4 part 1 for LTDC
HMD_function(eval = 'non-smooth')


##### MultiGait #####

MultiGait_function <- function(eval='smooth',python_df=TRUE){
  setwd("./gait")

  gait_df <- read.csv("gait.csv")

  gait_transformed <- gait_df %>%
    mutate(
      angle_rad = angle * (pi / 180),
      vec_x = cos(angle_rad),
      vec_y = sin(angle_rad)
    )
  
  list_of_lists <- gait_transformed %>%
    # UPDATED: group_by() must precede group_split()
    group_by(subject, condition, replication) %>%
    group_split() %>%
    lapply(function(trial_df) {
      
      inner_list <- trial_df %>%
        # UPDATED: group_by() must precede group_split()
        group_by(leg, joint) %>%
        group_split() %>%
        lapply(function(leg_joint_df) {
          
          leg_joint_df <- leg_joint_df[order(leg_joint_df$time), ]
          mat <- as.matrix(leg_joint_df[, c("vec_x", "vec_y")])
          rownames(mat) <- leg_joint_df$time
          
          return(mat)
        })
      
      # UPDATED: group_by() must precede group_keys()
      inner_keys <- trial_df %>% 
        group_by(leg, joint) %>%
        group_keys()
      
      names(inner_list) <- paste0(
        "Leg", inner_keys$leg, 
        "_Joint", inner_keys$joint
      )
      
      return(inner_list)
    })
  
  # UPDATED: group_by() must precede group_keys()
  outer_keys <- gait_transformed %>% 
    group_by(subject, condition, replication) %>%
    group_keys()
  
  names(list_of_lists) <- paste0(
    "Sub", outer_keys$subject, 
    "_Cond", outer_keys$condition, 
    "_Rep", outer_keys$replication
  )
  
  
  #####################################
  ## dissimilarity matrix
  list_of_matrices <- list_of_lists
  full_list <- list_of_matrices
  dissimilarity_matrix <- matrix(0,length(full_list),length(full_list))
  sub_name <- names(list_of_matrices)
  row.names(dissimilarity_matrix) = colnames(dissimilarity_matrix) = sub_name
  
  naive_dissimilarity_matrix <- dissimilarity_matrix
  SPD_dissimilarity_matrix <- dissimilarity_matrix  ##Spherical Procrustes distance
  ESD_dissimilarity_matrix <- dissimilarity_matrix ##Euclidean shape distance
  ED_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean distance
  EPD_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean Procrustes distance
  ##################################################
  
  if(eval=='smooth'){
    smooth_list <- list()
    
    for (i in 1:(length(full_list))) {
      sub_matrix_i <- full_list[[i]]
      sub_smooth_list <- list()
      for (k in 1:length(sub_matrix_i)) {
        sample_size <- nrow(sub_matrix_i[[k]])
        x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1)
        
        yin <- as.matrix(sub_matrix_i[[k]])
        row.names(yin) <- c(1:sample_size)
        trend_component <- LocSpheReg(xin = as.vector(x_input), yin = yin)
        bw_current <- trend_component$optns$bw 
        print(bw_current)
        sub_matrix_smooth <- trend_component$yout
        row.names(sub_matrix_smooth) <- row.names(sub_matrix_i[[k]])
        sub_smooth_list[[k]] <- sub_matrix_smooth
      }
      smooth_list[[i]] <- sub_smooth_list
      print(i)
    }
    names(smooth_list) <- sub_name
    
  } else{
    ## no smooth
    smooth_list <- list()
    
    for (i in 1:(length(full_list))) {
      sub_matrix_i <- full_list[[i]]
      smooth_list[[i]] <- sub_matrix_i
      print(i)
    }
    names(smooth_list) <- sub_name
  }
  
  
  ##################################################
  
  for (i in 1:(length(full_list)-1)) {
    for (j in (i+1):length(full_list)) {
      
      
      ##################################################
      ## smooth
      sub_matrix_i <- smooth_list[[i]]
      sub_matrix_j <- smooth_list[[j]]
      
      ##################################################
      ##shape
      dis_measure <- product_shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
      dissimilarity_matrix[i,j] <- dis_measure
      dissimilarity_matrix[j,i] <- dis_measure
      ##naive great circle distance
      dis_measure <- product_sphere_ts_dissimilarity_function(sub_matrix_i,sub_matrix_j)
      naive_dissimilarity_matrix[i,j] <- dis_measure
      naive_dissimilarity_matrix[j,i] <- dis_measure
      ##sphere Procrustes shape distance
      dis_measure <- D_rotProduct(sub_matrix_i,sub_matrix_j)
      SPD_dissimilarity_matrix[i,j] <- dis_measure
      SPD_dissimilarity_matrix[j,i] <- dis_measure
      ##Euclidean shape distance
      dis_measure <- product_function_shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
      ESD_dissimilarity_matrix[i,j] <- dis_measure
      ESD_dissimilarity_matrix[j,i] <- dis_measure
      
      ##Euclidean distance
      dis_measure <- product_function_dissimilarity_function(sub_matrix_i,sub_matrix_j)
      ED_dissimilarity_matrix[i,j] <- dis_measure
      ED_dissimilarity_matrix[j,i] <- dis_measure
      ##Euclidean Procrustes distance
      
      sub_matrix_i_save <- NULL
      sub_matrix_j_save <- NULL
      for (l in 1:length(sub_matrix_i)) {
        sub_matrix_i_save <- cbind(sub_matrix_i_save,sub_matrix_i[[l]])
        sub_matrix_j_save <- cbind(sub_matrix_j_save,sub_matrix_j[[l]])
      }
      dis_measure <- procrustes_distance(sub_matrix_i_save,sub_matrix_j_save)
      EPD_dissimilarity_matrix[i,j] <- dis_measure
      EPD_dissimilarity_matrix[j,i] <- dis_measure
    }
    print(i)
  }
  
  
  
  ##### classification
  
  all_items <- rownames(dissimilarity_matrix)
  
  all_labels <- sapply(strsplit(all_items, "_"), function(x) x[2])
  names(all_labels) <- all_items 
  
  set.seed(123) 
  
  n_total <- length(all_items)
  train_size <- floor(0.75 * n_total)
  
  # Randomly select items for the training set
  train_items <- sample(all_items, size = train_size)
  
  # The test items are whatever is left over
  test_items <- setdiff(all_items, train_items)
  
  # Subset the labels for the training set
  train_labels <- all_labels[train_items]
  
  # Keep the true labels for the test set so we can check our accuracy later
  true_test_labels <- all_labels[test_items]
  # ---------------------------------------------------------
  # 3. Define and Run the kNN Classifier for Distance Matrices
  # ---------------------------------------------------------
  result_matrix <- matrix(0,nrow = 2,ncol = 6)
  ##SSD
  results <- classification_evaluation_function(dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[1,1] <- round(results[[2]] * 100, 2)
  result_matrix[2,1] <- round(results[[4]] * 100, 2)
  
  
  ##SD
  results <- classification_evaluation_function(naive_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[1,2] <- round(results[[2]] * 100, 2)
  result_matrix[2,2] <- round(results[[4]] * 100, 2)
  
  ##SPD
  results <- classification_evaluation_function(SPD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[1,3] <- round(results[[2]] * 100, 2)
  result_matrix[2,3] <- round(results[[4]] * 100, 2)
  
  ##ESD
  results <- classification_evaluation_function(ESD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[1,4] <- round(results[[2]] * 100, 2)
  result_matrix[2,4] <- round(results[[4]] * 100, 2)
  
  
  ##ED
  results <- classification_evaluation_function(ED_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[1,5] <- round(results[[2]] * 100, 2)
  result_matrix[2,5] <- round(results[[4]] * 100, 2)
  
  ##ED
  results <- classification_evaluation_function(EPD_dissimilarity_matrix,
                                                all_items,
                                                train_items,
                                                all_labels,
                                                train_labels,
                                                test_items,
                                                true_test_labels) 
  
  result_matrix[1,6] <- round(results[[2]] * 100, 2)
  result_matrix[2,6] <- round(results[[4]] * 100, 2)
  
  
  ## prepare for Python 
  if (python_df==TRUE){
    names(smooth_list) <- country_labels$oecd_status[match(names(smooth_list), country_labels$ctr_name)]
    save(smooth_list,file = 'matrix_input.RData')

    dir.create("matrix_export", showWarnings = FALSE)

    for (i in seq_along(smooth_list)) {
      label_name <- names(smooth_list)[i]

      file_name <- paste0("matrix_export/", label_name, "_", i, ".csv")

      write.table(smooth_list[[i]], file = file_name, sep = ",", 
                  row.names = FALSE, col.names = FALSE)
    } 
  }
  return(result_matrix)
}

##Table 1 part 1 for MultiGait
MultiGait_function(eval = 'smooth')

##Table 4 part 1 for MultiGait
MultiGait_function(eval = 'non-smooth')
