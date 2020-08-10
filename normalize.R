# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# normalize the feature intensity values  

rm(list=ls())
library(tidyverse)
source("./lib/input.R")

help_text <- "
 NAME
    cluster.R

 SYNOPSIS
    cluster.R --input=<path_to_project_data>

 DESCRIPTION
    Normalize the intensity values of rt-mz features based on a random forest
    regression model.

 COMMAND LINE

    --input <path_to_project> [optional] (./data)
    
    --regex <file_pattern> [optional] (\\.ms1.fea.rds)

    --mdl <alignment_models_path> [optional] (NULL)

    --mis <min_missing> [optional 0-1] (0.75)

    --cpu <cpu_cores> [optional] (1)

 EXAMPLE

    Rscript cluster.R --input=./data --lctol=180

"

###############################################################################
# USER INPUT
ui_read_path        <- "./data"
ui_input_regex      <- "\\.ms1.fea.rds"
ui_save_models 		<- NULL
min_miss_ratio		<- 0.75
cpu_cores			<- 1

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--input", arg) ) ui_read_path <- arg_value
    if( grepl("--regex", arg) ) ui_input_regex <- arg_value
    if( grepl("--mdl", arg) ) ui_save_models <- arg_value
    if( grepl("--mis", arg) ) min_miss_ratio <- arg_value
    if( grepl("--cpu", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL

if(!is.null(message)) stop("ERROR\n", message)

file_list <- list.files(ui_read_path, 
                        pattern=ui_input_regex, 
                        recursive=T, 
                        full.names=T)

############################################################
# set up the parallization of FR ranger
set.seed(7224)
library(ranger)
library(doParallel)
cores 				<- cpu_cores
# getDoParWorkers() # Just checking
############################################################

d_big <- bigData(file_list)

# compute the means for lc mz, this will be our alignment seed
n_files <- df$file_name %>% unique() %>% length()
fat <- df %>%
    group_by(cluster_id) %>%
    summarise(
        n = length(cluster_id),
        mz = mean(mass_charge),
        lc = mean(elution_sec),
        int_log2 = mean(log2(intensity)),
        mz_ppm = sd(mass_charge)/mz * 1e6,
        lc_sd = sd(elution_sec),
        int_log2_rsd = sd(log2(intensity))/int_log2
    ) %>% filter(n >= n_files * min_miss_ratio)

mdf <- df %>% inner_join(fat, by="cluster_id") %>% 
    mutate(int_log2_mean = int_log2,
           int_log2_this = log2(intensity))

cat("normalize intensity values ...")
fit_list <- list()
for( i in 1:length(file_list) ){
    
    pb$tick()
    
    this_file_path <- file_list[i]
    this_file <- basename(this_file_path)
    
    stime <- Sys.time()
    
    dat_seed <- mdf %>% filter(file_name == this_file_name) 
    
    # loess regression
    # int_fit <- loess(int_log2_mean ~ int_log2_this, dat_seed, span=0.1)
    
    # random forest regression
    # int_fit <- randomForest(int_log2_mean ~ int_log2_this, dat_seed)
    
    cl <- makeCluster(cores)
    registerDoParallel(cores)
    int_fit <- ranger(int_log2_mean ~ int_log2_this, data = dat_seed)
    stopCluster(cl)
    
    dat_adjs <- file %>% readRDS() %>%
        mutate(int_log2_this = log2(intensity))
    
    int_pred <- predict(int_fit, data = dat_adjs)
    
    dat_adjs <- dat_adjs %>%
        mutate(intensity_norm = 2^int_pred$predictions) %>%
        mutate(intensity_norm = ifelse(is.na(intensity_norm), intensity, intensity_norm)) %>%
        select(-int_log2_this) %>%
        saveRDS(file = this_file_path)
    
    cat("\tdone", round(Sys.time() - stime,2), "sec\n")
    fit_list[[this_file_name]] <- list(int_fit = int_fit)  
}

if( !is.null(ui_save_models) )
    saveRDS(fit_list, file=ui_save_models)

