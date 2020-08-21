# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# read and merge data

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

readData <- function(file_name){
  if(grepl(".rds$", file_name)){
    df <- readRDS(file_name)
  } else if(grepl(".csv$", file_name)){
    df <- read.csv(file_name)
  }
    return(df)
}

saveData <- function(df, file_name){
  if(grepl(".rds$", file_name)){
    saveRDS(df, file_name)
  } else if(grepl(".csv$", file_name)){
    write.csv(df, file_name)
  }
}