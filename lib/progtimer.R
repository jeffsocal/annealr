# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# progress timer wrapper

library(progress)

progtimer <- function(size, text = "running ..."){
    cat(text, "\n")
    
    pb <- progress_bar$new(
        format = paste(size, ":percent :elapsed [:bar] :eta"),
        total = size, clear = FALSE, width= 60)
    
    return(pb)
}

hdf <- function(df){
    df %>% head() %>% as.data.frame()
}