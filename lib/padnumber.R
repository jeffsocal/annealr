# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# pad number helper

suppressMessages(
    library(
        tidyverse 
    )
)

num_pad <- Vectorize(
  function(x, w = NULL){
    if(is.null(w))
      w <- ceiling(log10(x))
    
    str_pad(x, w, "left", 0)
  }
)
