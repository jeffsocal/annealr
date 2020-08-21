# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster helper functions

suppressMessages(library(tidyverse))

# name all features for a given dataset
name_feaures <- function(df){
    df_n <- ceiling(log10(df %>% nrow()))
    df <- df %>%
        mutate(feature_id = paste0("F", str_pad(row_number(), df_n, "left","0")))
}

# create a global cluster reference
create_global_cluster <- function(df){
    df <- df %>% name_feaures() %>%
        mutate(cdist = 0) %>%
        select(elution_sec, mass_charge, charge, intensity, feature_id, cdist) %>%
        mutate(cluster_id = paste0("C", str_pad(row_number(), 6, "left","0")))
}