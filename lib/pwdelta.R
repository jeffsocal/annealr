# Jeff Jones
# SoCal Bioinformatics Inc. 2019
# 
# effecient all pairwise delta

suppressMessages(library(tidyverse))

pwdelta <- function(x_df=c(),
                    y_df=c(),
                    col_values=NULL,
                    col_names=NULL){
    
    x_df <- as.tibble(x_df) 
    y_df <- as.tibble(y_df) 
    
    x_vals <- x_df %>%
        dplyr::select_(col_values) %>%
        unlist()
    
    y_vals <- y_df %>%
        dplyr::select_(col_values) %>%
        unlist() 
    
    x_names <- x_df %>%
        dplyr::select_(col_names) %>%
        unlist() 
    
    y_names <- y_df %>%
        dplyr::select_(col_names) %>%
        unlist()
    
    x_vals_n <- x_df %>% nrow()
    y_vals_n <- y_df %>% nrow()
    
    if(x_vals_n != length(unique(x_names)))
        stop("Error: non-unique names for x")
    
    if(y_vals_n != length(unique(y_names)))
        stop("Error: non-unique names for y")
    
    xy_matrix <- array(x_vals, c(x_vals_n, y_vals_n))
    
    # matrix form
    xy_dif <- (t(xy_matrix) - y_vals) %>% as.numeric()
    
    x_row <- array(x_names, c(x_vals_n, y_vals_n)) %>% t() %>% as.character()
    y_row <- array(y_names, c(y_vals_n, x_vals_n)) %>% as.character()
    
    df <- tibble(
        ref=t(xy_matrix) %>% as.numeric(),
        dif=xy_dif,
        cluster_id=x_row,
        feature_id=y_row)
    
    return(df)
    
}