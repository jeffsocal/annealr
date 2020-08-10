# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster features across files  

rm(list=ls())
library(tidyverse)
options(warn=-1)
source("./lib/plotting.R")
source("./lib/clustering.R")
source('./lib/pairwise_delta.R')
source('./lib/progtimer.R')

help_text <- "
 NAME
    cluster.R

 SYNOPSIS
    cluster.R --input=<path_to_project_data>

 DESCRIPTION
    Clusters rt-mz features based on an optimized matrix evaluatiuon. Adds 
    the columns feature_id, cluster_id, and cdist which is the eluclidian 
    distance from the seed feature 

 COMMAND LINE

    --input <path_to_project> [optional] (./data)
    
    --regex <file_pattern> [optional] (\\.ms1.fea.rds)

    --mztol <mz_tolerance_in_daltons> [optional] (0.05)

    --lctol <lc_tolerance_in_sec> [optional] (30)

    --ztol <charge_state_tolerance> [optional] (0)

    --cpu <multi_threaded> [optional] (1)

    --chunk <matrix size> [optional] (64)

 EXAMPLE

    Rscript cluster.R --input=./data --lctol=180

"

###############################################################################
# USER INPUT
ui_read_path        <- "./data"
ui_input_regex      <- "\\.ms1.fea.rds"
ui_path_seed        <- NULL
mz_tol              <- 0.05          	# in daltons
lc_tol  		    <- 180           	# in seconds
cs_tol  		    <- 0            	# no tolerance on charge state
cpu_cores		    <- 1
chunk_size		    <- 64

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--input", arg) ) ui_read_path <- arg_value
    if( grepl("--regex", arg) ) ui_input_regex <- arg_value
    if( grepl("--mztol", arg) ) mz_tol <- arg_value
    if( grepl("--lctol", arg) ) lc_tol <- arg_value
    if( grepl("--ztol", arg) ) cs_tol <- arg_value
    if( grepl("--cpu", arg) ) cpu_cores <- arg_value
    if( grepl("--chunk", arg) ) chunk_size <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL

if(!is.null(message)) stop("ERROR\n", message)

cat("clustering parameters\n")
cat(" tolerance mz:                    ", mz_tol, "\n")
cat(" tolerance lc:                    ", lc_tol, "\n")
cat(" tolerance z:                     ", mz_tol, "\n")
cat(" cpus:                            ", cpu_cores, "\n")
cat(" chunk:                           ", chunk_size, "\n")

if( !is.null(ui_path_seed) ){
    d_gc <- ui_path_seed %>% readRDS()
    
    cat(" using a cluster seed file \n", ui_path_seed, "\n")   
}


file_list <- list.files(ui_read_path, 
                        pattern=ui_input_regex, 
                        recursive=T, 
                        full.names=T)

pb <- progtimer(length(file_list), "clustering ...")

for( i in 1:length(file_list) ){
    
    pb$tick()
    
    this_file <- file_list[i]
    
    df <- this_file %>% 
        readRDS() %>% 
        arrange(mass_charge) 
    
    #
    # remove previous clustering values
    #
    df_rm <- which(!grepl("cluster_id|cdist|feature_id", names(df)))
    df <- df[,df_rm]
    
    df_nr <- nrow(df)
    
    #
    # if no seed file exists create one
    #
    if( i == 1 && is.null(ui_path_seed) ){
        d_gc <- df %>% create_global_cluster()
        
        df <- df %>% full_join(
            d_gc, by=c('elution_sec', 'mass_charge', 'charge', 'intensity'))
        
        saveRDS(df, this_file)
        next()
    }
    
    ls_cf <- lcmz_cluster(
        tables = list('local' = df, 'global' = d_gc),
        mz_tol = mz_tol,
        lc_tol = lc_tol,
        cs_tol = cs_tol,
        cores = cpu_cores,
        chunk_size = chunk_size
    )
    
    d_gc <-ls_cf$global
    
    if(  df_nr != nrow(ls_cf$local))
        stop("\nwtf ... feature counts do not match!")
    
    saveRDS(ls_cf$local, this_file)
}
