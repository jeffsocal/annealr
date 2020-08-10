# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster helper functions

library(tidyverse)
library(parallel)

# cartiesian distance
cdist <- Vectorize(function(a, b){
  sqrt(a^2+b^2)
})

# name all features for a given dataset
name_feaures <- function(df){
  df_n <- ceiling(log10(df %>% nrow()))
  df <- df %>%
    mutate(feature_id = paste0("F", str_pad(row_number(), df_n, "left","0")))
}


create_global_cluster <- function(df){
  df <- df %>% name_feaures() %>%
    mutate(cdist = 0) %>%
    select(elution_sec, mass_charge, charge, intensity, feature_id, cdist) %>%
    mutate(cluster_id = paste0("C", str_pad(row_number(), 6, "left","0")))
}

# parallelizable function
lcmz_delta <- function(i_split, 
                       d_gc,
                       d_lc,
                       mz_tol  = 0.1,          # in daltons
                       lc_tol  = 60,           # in seconds
                       cs_tol  = 0            # no tolerance on charge state
){
  
  t_df <- d_lc[i_split,]
  t_d_gc <- d_gc %>%
    filter(mass_charge >= min(t_df$mass_charge) -1 &
             mass_charge <= max(t_df$mass_charge) +1 )
  
  # cluster by mass
  d_m <- pwdelta(t_d_gc, t_df, 'mass_charge', 'cluster_id') %>%
    filter(abs(dif) <= mz_tol)
  # cluster by elution
  d_e <- pwdelta(t_d_gc, t_df, 'elution_sec', 'cluster_id') %>%
    filter(abs(dif) <= lc_tol)
  # cluster by charge
  d_z <- pwdelta(t_d_gc, t_df, 'charge', 'cluster_id') %>%
    filter(dif == cs_tol) %>%
    inner_join(d_m, by=c("cluster_id", "feature_id"), suffix=c("", "_mass")) %>%
    inner_join(d_e, by=c("cluster_id", "feature_id"), suffix=c("", "_elution"))
  
  return(d_z)
}

lcmz_cluster <- function(
  # object containing the global and local clusters
  tables, 
  # clustering tolerances
  mz_tol  = 0.1,        # in daltons
  lc_tol  = 60,            # in seconds
  cs_tol  = 0,            # no tolerance on charge state
  chunk_size = 512,
  cores = 1
){
  
  d_lc <- tables$local
  d_gc <- tables$global
  
  # file to cluster with main-cluster-file
  d_lc <- d_lc %>% name_feaures() %>%
    mutate(cluster_id = feature_id)
  
  # divide the local finle into smaller chunks
  chunk_num = ceiling(nrow(d_lc)/chunk_size)
  v_split <- split(1:nrow(d_lc), sort(1:nrow(d_lc)%%chunk_num))
  
  n_df <- mclapply(v_split,
                   lcmz_delta,
                   d_gc, 
                   d_lc,
                   mz_tol  = mz_tol,          # in daltons
                   lc_tol  = lc_tol,           # in seconds
                   cs_tol  = cs_tol,            # no tolerance on charge state
                   mc.cores = cores) %>%
    bind_rows() %>% as_tibble() %>%
    select(-dif) %>%
    # compute the euclidean distance
    mutate(cdist = cdist(dif_mass, dif_elution)) %>%
    select(-c("dif_mass", "dif_elution"))
  
  d_clust <- d_lc %>%
    select(-cluster_id) %>%
    inner_join(
      n_df %>%
        # take closest cluster to feature center
        group_by(cluster_id) %>%
        top_n(1, -cdist) %>%
        ungroup() %>%
        # take closest feature to cluster center
        group_by(feature_id) %>%
        top_n(1, -cdist) %>%
        ungroup(),
      by="feature_id"
    )
  
  # get the max CF number
  max_cfn <- max(as.numeric(sub("C", "", d_gc$cluster_id)))
  
  # form a table of not-clustered
  d_clustnot <- d_lc %>%
    filter(!feature_id %in% d_clust$feature_id) %>%
    mutate(cluster_id = paste0("C", str_pad(row_number() + max_cfn, 6, "left","0"))) %>%
    mutate()
  
  # bind the not-clustered to the main-cluster-file
  d_gc <- d_gc %>%
    bind_rows(
      d_clustnot %>%
        mutate(cdist = 0) %>%
        select(elution_sec, mass_charge, charge, intensity, cluster_id, feature_id, cdist)
    )
  
  # final table of features all linked to cluster-references
  d_clust <- d_clust %>%
    bind_rows(
      d_clustnot
    )
  
  out <- list(
    'local' = d_clust,
    'global' = d_gc
  )
  
  return(out)
}
