# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# align clustered feature data  

rm(list=ls())
library(parallel)
suppressMessages(library(tidyverse))

help_text <- "
 NAME
    cluster.R

 SYNOPSIS
    cluster.R --input=<path_to_project_data>

 DESCRIPTION
    Align rt-mz features based on a loess regression fit model. 

 COMMAND LINE

    --input <path_to_project> [optional] (./data)
    
    --list <target_file_names> [optional] (NULL)
    
    --regex <file_pattern> [optional] (\\.ms1.fea.rds)

    --mdl <alignment_models_path> [optional] (NULL)

    --mis <min_missing> [optional 0-1] (0.75)

 EXAMPLE

    Rscript cluster.R --input=./data --lctol=180

"

###############################################################################
# USER INPUT
ui_read_path        <- "./data"
ui_input_regex      <- "\\.ms1.fea.rds"
ui_save_models 		<- NULL
ui_path_list        <- NULL
min_miss_ratio		<- 0.75
exe_path            <- "./"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--input", arg) ) ui_read_path <- arg_value
    if( grepl("--list", arg) ) ui_path_list <- arg_value
    if( grepl("--regex", arg) ) ui_input_regex <- arg_value
    if( grepl("--mdl", arg) ) ui_save_models <- arg_value
    if( grepl("--mis", arg) ) min_miss_ratio <- arg_value
    if( grepl("--exp", arg) ) exe_path <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

min_miss_ratio <- 1 - min_miss_ratio

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


if( !is.null(ui_path_list) ){
    file_select <- ui_path_list %>% read.csv()
    
    w <- which(grepl(
        paste(file_select[,1], collapse = "|"), 
        file_list))    
    if( length(w) > 0 )
        file_list <- file_list[w]
}

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
        mz_ppm = sd(mass_charge)/mz * 1e6,
        lc_sd = sd(elution_sec),
        it_rsd = sd(intensity)/it,
        .groups = 'drop'
    ) %>%
    filter(n >= n_files * min_miss_ratio) 

if( nrow(fat) < 100 ) {
    stop(" too few clusters remain after filtering for ", min_miss_ratio * 100, "% completeness")    
}

# merge the means back into the data
mdf <- d_big %>% inner_join(fat, by="cluster_id")

fit_list <- list()

cat(" aligning individual scans ... \n")
n_files <- length(file_list)
for( i in 1:n_files ){
    
    this_file_path <- file_list[i]
    this_file <- basename(this_file_path)
    
    cat("  ", i, "of", n_files, this_file, "\n")
    
    # define the seed  
    dat_seed <- mdf %>% filter(file_name == this_file)
    
    # non-linear regression loess
    mz_fit <- loess(mz ~ mass_charge, dat_seed, span=0.05, normalize = F)
    lc_fit <- loess(lc ~ elution_sec, dat_seed, span=0.05, normalize = F)
    
    # non-linear regression spline
    # mz_fit <- lm(mz ~ bs(mass_charge, df=5), dat_seed, normalize = F)
    # lc_fit <- lm(lc ~ bs(elution_sec, df=15), dat_seed, normalize = F)
    
    dat_adjs <- this_file_path %>% readRDS() %>%
        select(-matches("elution_sec_org|mass_charge_org")) %>%
        mutate(mass_charge_fit = predict(mz_fit, newdata = .)) %>%
        mutate(elution_sec_fit = predict(lc_fit, newdata = .)) %>%
        mutate(mass_charge_fit = ifelse(is.na(mass_charge_fit), mass_charge, mass_charge_fit)) %>%
        mutate(elution_sec_fit = ifelse(is.na(elution_sec_fit), elution_sec, elution_sec_fit)) %>%
        rename(
            mass_charge_org = mass_charge,
            elution_sec_org = elution_sec,
            mass_charge = mass_charge_fit,
            elution_sec = elution_sec_fit
        ) %>% 
        saveRDS(file = this_file_path)
    
    fit_list[[this_file]] <- list(mz_fit = mz_fit, lc_fit = lc_fit)  
}

if( !is.null(ui_save_models) )
    saveRDS(fit_list, file=ui_save_models)

cat("completed\n")
