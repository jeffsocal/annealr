# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# normalize the feature intensity values  

rm(list=ls())
suppressMessages(library(tidyverse))

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
min_miss_ratio		<- 0.50
cpu_cores			<- 1
exe_path            <- "./"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--input", arg) ) ui_read_path <- arg_value
    if( grepl("--regex", arg) ) ui_input_regex <- arg_value
    if( grepl("--mdl", arg) ) ui_save_models <- arg_value
    if( grepl("--mis", arg) ) min_miss_ratio <- arg_value
    if( grepl("--exp", arg) ) exe_path <- arg_value
    if( grepl("--cpu", arg) ) cpu_cores <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

source(paste0(exe_path, "/lib/input.R"))
source(paste0(exe_path, "/lib/progtimer.R"))
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
# set.seed(7224)
# library(ranger)
# library(doParallel)
# cores 				<- cpu_cores
# getDoParWorkers() # Just checking
############################################################

n_files <- file_list %>% length()
d_big <- bigData(file_list, " reading and combining all data sets ...")

# compute the means for lc mz, this will be our alignment seed
fat <- d_big %>%
    group_by(cluster_id) %>%
    summarise(
        n = length(cluster_id),
        mz = mean(mass_charge),
        lc = mean(elution_sec),
        it = mean(intensity),
        it_sd = sd(intensity),
        int_log2 = mean(log2(intensity)),
        mz_ppm = sd(mass_charge)/mz * 1e6,
        lc_sd = sd(elution_sec),
        int_log2_rsd = sd(log2(intensity))/int_log2,
        .groups = "drop"
    ) %>% 
    filter(n >= n_files * min_miss_ratio)

mdf <- d_big %>% 
    inner_join(fat, by="cluster_id") %>% 
    mutate(int_log2_mean = int_log2,
           int_log2_this = log2(intensity))

fit_list <- list()

cat(" normalizing intensity values ... \n")
n_files <- length(file_list)
for( i in 1:n_files ){
    
    this_file_path <- file_list[i]
    this_file <- basename(this_file_path)
    
    cat("  ", i, "of", n_files, this_file, "...")
    
    stime <- Sys.time()
    
    dat_seed <- mdf %>% filter(file_name == this_file) 
    
    # loess regression
    int_fit <- loess(int_log2_mean ~ int_log2_this, dat_seed, span=0.1)
    
    # random forest regression
    # int_fit <- randomForest(int_log2_mean ~ int_log2_this, dat_seed)
    
    # ranger package
    # cl <- makeCluster(cores)
    # registerDoParallel(cores)
    # int_fit <- ranger(int_log2_mean ~ int_log2_this, data = dat_seed)
    # stopCluster(cl)
    
    dat_adjs <- this_file_path %>% 
        readRDS() %>%
        mutate(int_log2_this = log2(intensity))
    
    int_pred <- predict(int_fit, newdata = dat_adjs)

    # ranger
    # int_pred <- int_pred$predictions
    
    dat_adjs <- dat_adjs %>%
        mutate(intensity_norm = 2^int_pred) %>%
        mutate(intensity_norm = ifelse(is.na(intensity_norm), intensity, intensity_norm)) %>%
        select(-int_log2_this) %>%
        saveRDS(file = this_file_path)
    
    cat("\tdone", round(Sys.time() - stime,2), "sec\n")
    fit_list[[this_file]] <- list(int_fit = int_fit)  
}

if( !is.null(ui_save_models) )
    saveRDS(fit_list, file=ui_save_models)

cat("assessment:")
n_files <- file_list %>% length()
d_big <- bigData(file_list, " reading and combining all data sets ...")
# compute the means for lc mz, this will be our alignment seed
fatp <- d_big %>%
    group_by(cluster_id) %>%
    summarise(
        n = length(cluster_id),
        mz = mean(mass_charge),
        lc = mean(elution_sec),
        it = mean(intensity_norm),
        mz_sd = sd(mass_charge),
        lc_sd = sd(elution_sec),
        it_sd = sd(intensity_norm),
        .groups = 'drop'
    )

fat <- fat %>% filter(n >= n_files * .50)
fatp <- fatp %>% filter(n >= n_files * .50)

cat(" Int accuracy 50% inc PRE\n")
cat("  CV mean              ", (fat$it_sd / fat$it) %>% mean(rm.na=T) %>% signif(5), "\n")
cat(" Int accuracy 50% inc POST\n")
cat("  CV mean              ", (fatp$it_sd / fatp$it) %>% mean(rm.na=T) %>% signif(5), "\n")
