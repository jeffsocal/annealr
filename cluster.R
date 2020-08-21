# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster features across files  

rm(list=ls())
suppressMessages(library(tidyverse))
options(warn=-1)

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

    --seedfile <path_to_seed> [optional] (NULL)
    
    --list <target_file_name> [optional] (NULL)
    
    --regex <file_pattern> [optional] (\\.ms1.fea.rds)

    --mztol <mz_tolerance_in_daltons> [optional] (0.05)

    --lctol <lc_tolerance_in_sec> [optional] (30)

    --ztol <charge_state_tolerance> [optional] (0)

    --cpu <multi_threaded> [optional] (1)

    --chunk <matrix size> [optional] (64)

    --exp <path to executable> [optional] (./)

 EXAMPLE

    Rscript cluster.R --input=./data --lctol=180

"

###############################################################################
# USER INPUT
ui_read_path        <- "./data"
ui_input_regex      <- "\\.ms1.fea.rds"
ui_path_seed        <- NULL
ui_path_list        <- NULL
mz_tol              <- 0.05          	# in daltons
lc_tol  		    <- 180           	# in seconds
cs_tol  		    <- 0            	# no tolerance on charge state
cpu_cores		    <- 1
chunk_size		    <- 64
exe_path            <- "."

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--input", arg) ) ui_read_path <- arg_value
    if( grepl("--regex", arg) ) ui_input_regex <- arg_value
    if( grepl("--list", arg) ) ui_path_list <- arg_value
    if( grepl("--seed", arg) ) ui_path_seed <- arg_value
    if( grepl("--mztol", arg) ) mz_tol <- arg_value
    if( grepl("--lctol", arg) ) lc_tol <- arg_value
    if( grepl("--ztol", arg) ) cs_tol <- arg_value
    if( grepl("--chunk", arg) ) chunk_size <- arg_value
    if( grepl("--exp", arg) ) exe_path <- arg_value
    if( grepl("--cpu", arg) ) cpu_cores <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

source(paste0(exe_path, "/lib/cdist.R"))
source(paste0(exe_path, "/lib/input.R"))
source(paste0(exe_path, "/lib/pwdelta.R"))
source(paste0(exe_path, "/lib/mzedelta.R"))
source(paste0(exe_path, "/lib/mzecluster.R"))
source(paste0(exe_path, "/lib/progtimer.R"))
source(paste0(exe_path, "/lib/references.R"))
###############################################################################
# INPUT VALIDATION
message <- NULL

if(!is.null(message)) stop("ERROR\n", message)

cat("clustering parameters\n")
cat(" tolerance mz:                    ", mz_tol, "\n")
cat(" tolerance lc:                    ", lc_tol, "\n")
cat(" tolerance z:                     ", cs_tol, "\n")
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

if( !is.null(ui_path_list) ){
    file_select <- ui_path_list %>% read.csv()
    
    w <- which(grepl(
        paste(file_select[,1], collapse = "|"), 
        file_list))    
    if( length(w) > 0 )
        file_list <- file_list[w]
}

cat("clustering ... \n")
n_files <- length(file_list)
for( i in 1:n_files ){
    
    this_file_path <- file_list[i]
    this_file <- basename(this_file_path)
    
    cat("  ", i, "of", n_files, this_file, "\n")
    df <- this_file_path %>% 
        readRDS() %>% 
        select(elution_sec, mass_charge, charge, intensity,
               quality_overall, quality_elution_sec, quality_mass_charge,
               FWHM, label, score_correlation, score_fit,
               spectrum_index, spectrum_native_id,
               feature_index) %>%
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
            d_gc, by=c('elution_sec', 'mass_charge', 'charge', 'intensity')) %>%
            select(-c("feature_id", "cdist"))
        
        saveRDS(df, this_file_path)
        next()
    }
    
    ls_cf <- mzecluster(
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
    
    saveRDS(
        ls_cf$local %>%
            select(-matches("feature_id|cdist|ref_*|dif_*"))
    , this_file_path)
}

cat("completed\n")

cat("assessment:")
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
        mz_sd = sd(mass_charge),
        lc_sd = sd(elution_sec),
        it_rsd = sd(intensity)/it,
        .groups = 'drop'
    ) 

cat(" n features \n")
cat("  100% inc             ", fat %>% filter(n >= n_files * 1) %>% nrow(), "\n")
cat("  75% inc              ", fat %>% filter(n >= n_files * .75) %>% nrow(), "\n")
cat("  50% inc              ", fat %>% filter(n >= n_files * .50) %>% nrow(), "\n")
cat("  25% inc              ", fat %>% filter(n >= n_files * .25) %>% nrow(), "\n")
cat("  Total                ", fat %>% nrow(), "\n")

fat <- fat %>% filter(n >= n_files * .50)
cat(" mz accuracy 50% inc \n")
cat("  95% CI               ", (fat$mz_sd %>% sd() * 2.56) %>% signif(3), "\n")
cat("  CV mean              ", (fat$mz_sd / fat$mz) %>% mean(rm.na=T) %>% signif(3), "\n")
cat(" rt accuracy 50% inc \n")
cat("  95% CI               ", (fat$lc_sd %>% sd() * 2.56) %>% signif(3), "\n")
cat("  CV mean              ", (fat$lc_sd / fat$lc) %>% mean(rm.na=T) %>% signif(3), "\n")