# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# read and merge data

source('./lib/progtimer.R')

bigData <- function(file_paths = NULL, message = "merging data"){
    df <- c()
    pb <- progtimer(length(file_paths), message)
    for( this_file_path in file_paths ){
        
        this_file <- basename(this_file_path)
        
        pb$tick()
        df <- df %>% 
            bind_rows(this_file_path %>% 
                          readRDS() %>%
                          mutate(file_name = this_file))
    }
    return(df)
}