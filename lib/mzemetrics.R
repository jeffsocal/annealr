# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# metrics calculations

suppressMessages(library(tidyverse))

metricChrom <- function(df){
  
}

metricMSn <- function(df, fraction = c('pep_z', 'file_name')){
  
  library(stringr)
  library(reshape2)
  
  # over 5 min windows
  #   s/n 
  #   msn quality
  
  df$pre_rt_5min <- round(df$pre_rt * 5)/5
  
  fraction <- intersect(fraction, c('pre_rt_5min',colnames(df)))
  
  idn <- ddply(df, fraction, summarize,
               msn_sn_median = median(msn_sn),
               msn_qual_median = median(msn_qual)
  )
  
  return(idn)  
  
}

metricIDs <- function(df, fraction = c('pep_z', 'file_name')){
  
  library(stringr)
  library(reshape2)
  
  # count of IDs by charge state
  #   total, unique, by length
  # summary of ID scores
  #   median, IQR
  # ratio of internal k|r unique / total unique
  
  df$pep_length <- str_length(df$pep_seq)
  
  fraction <- intersect(fraction, colnames(df))
  
  idn <- ddply(df, fraction, summarize,
               pep_total_n = length(pep_seq),
               pep_unique_n = length(unique(pep_seq)),
               pep_mcleave_n = length(unique(pep_seq[grepl("R|K", pep_seq)])),
               pep_len_median = median(pep_length),
               pep_len_qnt20 = qunatile(pep_length,c(.2))[1],
               pep_len_qnt80 = qunatile(pep_length,c(.8))[1],
               id_log2pval_median = median(id_pvalue),
               id_log2pval_iqr = IQR(id_pvalue)
  )
  
  return(idn)  
  
}