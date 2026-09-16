library(movMF)
library(trust)
library(parallel)
library(foreach)
library(doSNOW)

source('basic_function.R')



sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  
  for (d_current in c(6,101)) {
    for (cluster_num in c(30,90)) {
      full_list <- shaple_cluster_sim(seed=rep,function_size=cluster_num,sample.size=200,d_set=d_current)
      #####################################
      ## dissimilarity matrix
      dissimilarity_matrix <- matrix(0,length(full_list),length(full_list))
      cluster_name <- names(full_list)
      row.names(dissimilarity_matrix) = colnames(dissimilarity_matrix) = 1:nrow(dissimilarity_matrix)
      
      naive_dissimilarity_matrix <- dissimilarity_matrix
      SPD_dissimilarity_matrix <- dissimilarity_matrix  ##Spherical Procrustes distance
      ESD_dissimilarity_matrix <- dissimilarity_matrix ##Euclidean shape distance
      ED_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean distance
      EPD_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean Procrustes distance
      
      ##################################################
      
      ## smooth
      smooth_list <- list()
      
      for (i in 1:(length(full_list))) {
        sub_matrix_i <- full_list[[i]]
        
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
      names(smooth_list) <- cluster_name
      
      
      ##################################################
      
      for (i in 1:(length(full_list)-1)) {
        for (j in (i+1):length(full_list)) {
          sub_matrix_i <- smooth_list[[i]]
          sub_matrix_j <- smooth_list[[j]]
          
          ##################################################
          ##shape
          dis_measure <- shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
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
      dissimilarity_matrix_list <- list(
        dissimilarity_matrix,
        naive_dissimilarity_matrix,
        SPD_dissimilarity_matrix,
        ESD_dissimilarity_matrix,
        ED_dissimilarity_matrix,
        EPD_dissimilarity_matrix
      )
      result_list[[list_ind]] <- dissimilarity_matrix_list
      list_ind <- list_ind + 1
    }
  }
  
  return(result_list)
}

no_smooth_sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  
  for (d_current in c(6,101)) {
    for (cluster_num in c(30,90)) {
      full_list <- shaple_cluster_sim(seed=rep,function_size=cluster_num,sample.size=200,d_set=d_current)
      #####################################
      ## dissimilarity matrix
      dissimilarity_matrix <- matrix(0,length(full_list),length(full_list))
      cluster_name <- names(full_list)
      row.names(dissimilarity_matrix) = colnames(dissimilarity_matrix) = 1:nrow(dissimilarity_matrix)
      
      naive_dissimilarity_matrix <- dissimilarity_matrix
      SPD_dissimilarity_matrix <- dissimilarity_matrix  ##Spherical Procrustes distance
      ESD_dissimilarity_matrix <- dissimilarity_matrix ##Euclidean shape distance
      ED_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean distance
      EPD_dissimilarity_matrix <- dissimilarity_matrix  ##Euclidean Procrustes distance
      
      
      for (i in 1:(length(full_list)-1)) {
        for (j in (i+1):length(full_list)) {
          
          ##################################################
          ## non-smooth
          sub_matrix_i <- full_list[[i]]
          sub_matrix_j <- full_list[[j]]

          ##shape
          dis_measure <- shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
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
        #print(i)
      }
      dissimilarity_matrix_list <- list(
        dissimilarity_matrix,
        naive_dissimilarity_matrix,
        SPD_dissimilarity_matrix,
        ESD_dissimilarity_matrix,
        ED_dissimilarity_matrix,
        EPD_dissimilarity_matrix
      )
      result_list[[list_ind]] <- dissimilarity_matrix_list
      list_ind <- list_ind + 1
    }
  }
  
  return(result_list)
}




# Setup Cluster
nCores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(nCores)
registerDoSNOW(cl)

iterations <- 200

# Progress Bar
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# Parallel Loop
sim_results <- foreach(
  i = 1:iterations,
  .combine = 'c',               
  .multicombine = TRUE,           
  .options.snow = opts,
  .packages = c("movMF","trust")          
) %dopar% {
  source('basic_function.R') 
  
  # Ensure the return is wrapped in a list so 'c' combines them into a list of 200
  ##list(sim_function(rep = i))
  list(no_smooth_sim_function(rep = i))
}

parallel::stopCluster(cl)
close(pb)


accuracy_matrix <- matrix(0,nrow = 12,ncol = 4)

for (j in 1:4) {
  for (i in 1:200){
    dissimilarity_matrix_list <- sim_results[[i]][[j]]
    dissimilarity_matrix <- dissimilarity_matrix_list[[1]]
    naive_dissimilarity_matrix <- dissimilarity_matrix_list[[2]]
    SPD_dissimilarity_matrix <- dissimilarity_matrix_list[[3]]
    ESD_dissimilarity_matrix <- dissimilarity_matrix_list[[4]]
    ED_dissimilarity_matrix <- dissimilarity_matrix_list[[5]]
    EPD_dissimilarity_matrix <- dissimilarity_matrix_list[[6]]

    all_items <- rownames(dissimilarity_matrix)
    all_labels <- rep(1:3, length.out = length(all_items))
    names(all_labels) <- all_items

    set.seed(123) 
    
    n_total <- length(all_items)
    train_size <- floor(0.75 * n_total)
    
    train_items <- sample(all_items, size = train_size)
    test_items <- setdiff(all_items, train_items)
    train_labels <- all_labels[train_items]
    true_test_labels <- all_labels[test_items]

    ##SSD
    results <- classification_evaluation_function(dissimilarity_matrix,
                                                  all_items,
                                                  train_items,
                                                  train_labels,
                                                  test_items,
                                                  true_test_labels) 
    
    accuracy_matrix[1,j] <- accuracy_matrix[1,j] + results[[2]] ## random
    accuracy_matrix[7,j] <- accuracy_matrix[5,j] + results[[4]] ## LOOCV
    
    
    ##SD
    results <- classification_evaluation_function(naive_dissimilarity_matrix,
                                                  all_items,
                                                  train_items,
                                                  train_labels,
                                                  test_items,
                                                  true_test_labels) 
    
    accuracy_matrix[2,j] <- accuracy_matrix[2,j] + results[[2]] ## random
    accuracy_matrix[8,j] <- accuracy_matrix[6,j] + results[[4]] ## LOOCV
    
    ##SPD
    results <- classification_evaluation_function(SPD_dissimilarity_matrix,
                                                  all_items,
                                                  train_items,
                                                  train_labels,
                                                  test_items,
                                                  true_test_labels) 
    
    accuracy_matrix[3,j] <- accuracy_matrix[2,j] + results[[2]] ## random
    accuracy_matrix[9,j] <- accuracy_matrix[6,j] + results[[4]] ## LOOCV
    
    ##ESD
    results <- classification_evaluation_function(ESD_dissimilarity_matrix,
                                                  all_items,
                                                  train_items,
                                                  train_labels,
                                                  test_items,
                                                  true_test_labels) 
    
    accuracy_matrix[4,j] <- accuracy_matrix[3,j] + results[[2]] ## random
    accuracy_matrix[10,j] <- accuracy_matrix[7,j] + results[[4]] ## LOOCV
    
    ##ED
    results <- classification_evaluation_function(ED_dissimilarity_matrix,
                                                  all_items,
                                                  train_items,
                                                  train_labels,
                                                  test_items,
                                                  true_test_labels) 
    
    accuracy_matrix[5,j] <- accuracy_matrix[4,j] + results[[2]] ## random
    accuracy_matrix[11,j] <- accuracy_matrix[8,j] + results[[4]] ## LOOCV
    
    ##EPD
    results <- classification_evaluation_function(EPD_dissimilarity_matrix,
                                                  all_items,
                                                  train_items,
                                                  train_labels,
                                                  test_items,
                                                  true_test_labels) 
    
    accuracy_matrix[6,j] <- accuracy_matrix[4,j] + results[[2]] ## random
    accuracy_matrix[12,j] <- accuracy_matrix[8,j] + results[[4]] ## LOOCV
    
    
  }
}



full_results <- round(t(accuracy_matrix) / 200 * 100, 2)
rbind(full_results[,1:6],full_results[,7:12])

####################################################################################################
####################################################################################################
## prepare for Python 
####################################################################################################
####################################################################################################

## smooth
smooth_list_function <- function(rep=1,cluster_num=30,d_current=6){
  full_list <- shaple_cluster_sim(seed=rep,function_size=cluster_num,sample.size=200,d_set=d_current)
  cluster_name <- names(full_list)
  ## smooth
  smooth_list <- list()
  
  for (i in 1:(length(full_list))) {
    sub_matrix_i <- full_list[[i]]
    
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
  names(smooth_list) <- cluster_name
  return(smooth_list)
}

no_smooth_list_function <- function(rep=1,cluster_num=30,d_current=6){
  full_list <- shaple_cluster_sim(seed=rep,function_size=cluster_num,sample.size=200,d_set=d_current)
  return(full_list)
}


##cluster_num=30,d_current=6
##cluster_num=90,d_current=6
##cluster_num=30,d_current=101
##cluster_num=90,d_current=101



par_list <- list(c(30,6),c(90,6),c(30,101),c(90,101))
dir_list <- c("sim_data/matrix_export1","sim_data/matrix_export2","sim_data/matrix_export3","sim_data/matrix_export4")
for (ind in 1:4) {
  smooth_list <- smooth_list_function(cluster_num=par_list[[ind]][1],d_current=par_list[[ind]][2])
  ##smooth_list <- no_smooth_list_function(cluster_num=par_list[[ind]][1],d_current=par_list[[ind]][2])
  # Create a folder to hold the data
  dir.create(dir_list[ind], showWarnings = FALSE)
  
  # Loop through the list and save each matrix
  for (i in seq_along(smooth_list)) {
    label_name <- names(smooth_list)[i]
    file_name <- paste0(dir_list[ind],"/", label_name, "_", i, ".csv")
    write.table(smooth_list[[i]], file = file_name, sep = ",", 
                row.names = FALSE, col.names = FALSE)
  }
  
}


#################################
#################################
## computation time
#################################
#################################
cluster_num <- 90
d_current <- 6
full_list <- shaple_cluster_sim(seed=1,function_size=cluster_num,sample.size=200,d_set=d_current)

##### smoothing#####
start.time <- Sys.time()
sub_matrix_i <- full_list[[1]]

sample_size <- nrow(sub_matrix_i)
x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1)

trend_component <- LocSpheReg(xin = as.vector(x_input), yin = sub_matrix_i)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
## m=6 27.8072 secs
## m=101 1.648414 mins


##### classification#####
start.time <- Sys.time()
dissimilarity_matrix <- matrix(0,length(full_list),length(full_list))
cluster_name <- names(full_list)
row.names(dissimilarity_matrix) = colnames(dissimilarity_matrix) = 1:nrow(dissimilarity_matrix)
for (i in 1:(length(full_list)-1)) {
  for (j in (i+1):length(full_list)) {
    sub_matrix_i <- full_list[[i]]
    sub_matrix_j <- full_list[[j]]
    dis_measure <- shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
    ##dis_measure <- sphere_ts_dissimilarity_function(sub_matrix_i,sub_matrix_j)
    ##dis_measure <- sphere_procrustes_distance(sub_matrix_i,sub_matrix_j)
    ##dis_measure <- multi_function_dissimilarity_function(sub_matrix_i,sub_matrix_j)
    ##dis_measure <- function_shape_dissimilarity_function_1D(sub_matrix_i,sub_matrix_j)
    ##dis_measure <- procrustes_distance(sub_matrix_i,sub_matrix_j)
    dissimilarity_matrix[i,j] <- dis_measure
    dissimilarity_matrix[j,i] <- dis_measure
  }
  print(i)
}
all_items <- rownames(dissimilarity_matrix)
all_labels <- rep(1:3, length.out = length(all_items))
names(all_labels) <- all_items
set.seed(123) 
n_total <- length(all_items)
train_size <- floor(0.75 * n_total)
train_items <- sample(all_items, size = train_size)
test_items <- setdiff(all_items, train_items)
train_labels <- all_labels[train_items]
true_test_labels <- all_labels[test_items]
results <- classification_evaluation_function(dissimilarity_matrix,
                                              all_items,
                                              train_items,
                                              train_labels,
                                              test_items,
                                              true_test_labels) 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
## PSSD
## n=30 m=6 1.765035 secs
## n=90 m=6  15.89417 secs
## n=30 m=101 26.10263 secs
## n=90 m=101 3.880579 mins

## SD
## n=30 m=6 0.04341984 secs
## n=90 m=6 0.106971 secs
## n=30 m=101 0.05094886 secs
## n=90 m=101 0.1650529 secs

## SPD
## n=30 m=6 2.70195 secs
## n=90 m=6 24.64281 secs
## n=30 m=101  1.25928 mins
## n=90 m=101  12.03516 mins

## ESD
## n=30 m=6 0.02816796 secs
## n=90 m=6 0.08910704 secs
## n=30 m=101 0.04875088 secs 
## n=90 m=101 0.2400329 secs

## ED
## n=30 m=6 0.04805613 secs
## n=90 m=6 0.1509609 secs
## n=30 m=101  0.124294 secs
## n=90 m=101  0.814467 secs

## EPD
## n=30 m=6 0.226953 secs
## n=90 m=6 1.432196 secs
## n=30 m=101  2.944537 secs
## n=90 m=101  26.68911 secs





