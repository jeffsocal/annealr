# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# impute the missing feature intensity values 

rm(list=ls())
suppressMessages(library(tidyverse))

help_text <- "
 NAME
    cluster.R

 SYNOPSIS
    cluster.R --input=<path_to_project_data>

 DESCRIPTION
    Impute missing intensity values for rt-mz features based on a random forest
    model.

 COMMAND LINE

    --input <path_to_project> [optional] (./data)
    
    --regex <file_pattern> [optional] (\\.ms1.fea.rds)

    --out <labelfree_matrix> [optional] (labelfree_matrix)

    --impute <TRIUE/FALSE> [optional] (TRUE)

    --cpu <cpu_cores> [optional] (1)

 EXAMPLE

    Rscript cluster.R --input=./data --lctol=180

"

###############################################################################
# USER INPUT
ui_read_path        <- "./data"
ui_input_regex      <- "\\.ms1.fea.rds"
ui_save_matrix      <- "labelfree_matrix"
ui_impute_yn        <- TRUE
cpu_cores			  <- 1
exe_path            <- "./"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--input", arg) ) ui_read_path <- arg_value
    if( grepl("--regex", arg) ) ui_input_regex <- arg_value
    if( grepl("--out", arg) ) ui_save_matrix <- arg_value
    if( grepl("--impute", arg) ) ui_impute_yn <- arg_value
    if( grepl("--exp", arg) ) exe_path <- arg_value
    if( grepl("--cpu", arg) ) cpu_cores <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

source(paste0(exe_path, "lib/input.R"))
source(paste0(exe_path, "lib/progtimer.R"))
###############################################################################
# INPUT VALIDATION
message <- NULL

if(!is.null(message)) stop("ERROR\n", message)

file_list <- list.files(ui_read_path, 
                        pattern=ui_input_regex, 
                        recursive=T, 
                        full.names=T)

d_big <- bigData(file_list)

# create the matrix -- rf-norm
d_mx <- d_big %>%
    select(file_name, cluster_id, intensity_norm) %>%
    spread(cluster_id, intensity_norm)

d_mx <- d_mx %>% as.data.frame()
rownames(d_mx) <- d_mx$file_name
m_can <- d_mx %>% select(-file_name) %>% as.matrix()

if(ui_impute_yn == TRUE) {
    
    # impute the missing values
    library(missForest)
    library(doParallel)
    set.seed(7224)
    cores <- cpu_cores
    
    cat("create clustering cores ... \n")
    cl <- makeCluster(cores)
    registerDoParallel(cores)
    
    cat("  impute rf matrix ... \n")
    # should be one of “no”, “variables”, “forests”
    m_ican <- missForest(m_can, parallelize = "variables")
    
    
    stopCluster(cl)
}

saveRDS(m_ican, paste0(ui_read_path, sub("\\.rds$", "", ui_save_matrix), ".rds"))
cat("DONE\n")
