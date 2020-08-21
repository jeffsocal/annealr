# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster helper functions

suppressMessages(library(tidyverse))

# parallelizable function
mzedelta <- function(i_split, 
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
        filter(abs(dif) <= cs_tol) %>%
        rename(dif_charge = dif, ref_charge = ref) %>%
        inner_join(d_m, by=c("cluster_id", "feature_id")) %>%
        rename(dif_mass = dif, ref_mass = ref) %>%
        inner_join(d_e, by=c("cluster_id", "feature_id")) %>%
        rename(dif_elution = dif, ref_elution = ref)
    
    return(d_z)
}
