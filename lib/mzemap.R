# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster helper functions

library(parallel)
suppressMessages(library(tidyverse))

mzemap <- function(
    # object containing the reference and multiples
    tables, 
    # clustering tolerances
    mz_tol  = 0.1,          # in daltons
    lc_tol  = 60,           # in seconds
    cs_tol  = 0,            # no tolerance on charge state
    chunk_size = 512,
    cores = 1
){
    
    d_gc <- tables$reference
    d_lc <- tables$multiples
    
    d_gc <- d_gc %>% 
        rename(cluster_id = feature_index)
    
    d_lc <- d_lc %>% 
        rename(feature_id = scan_index) %>%
        mutate(cluster_id = feature_id,
               elution_sec = rt_time_sec,
               mass_charge = precursor_mz,
               charge = precursor_z,
               intensity = precursor_int)
    
    # divide the local file into smaller chunks
    chunk_num = ceiling(nrow(d_lc)/chunk_size)
    v_split <- split(1:nrow(d_lc), sort(1:nrow(d_lc)%%chunk_num))
    
    n_df <- mclapply(v_split,
                     mzedelta,
                     d_gc, 
                     d_lc,
                     mz_tol  = mz_tol,          # in daltons
                     lc_tol  = lc_tol,           # in seconds
                     cs_tol  = cs_tol,            # no tolerance on charge state
                     mc.cores = cores) %>%
        bind_rows() %>% as_tibble() %>%
        mutate(dif_mass_ppm = dif_mass / ref_mass * 1e6) %>%
        # compute the euclidean distance
        mutate(cdist = cdist(dif_mass_ppm, dif_elution))
    
    d_clust <- d_lc %>%
        select(-cluster_id) %>%
        inner_join(
            n_df %>%
                # take closest reference feature 
                # to center of precursor mz/rt
                group_by(feature_id) %>%
                top_n(1, -cdist) %>%
                ungroup(),
            by="feature_id"
        ) %>% 
        select(feature_index = cluster_id,
               scan_index = feature_id,
               dif_mass, dif_elution)
    
    return(d_clust)
}