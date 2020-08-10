# Jeff Jones
# SoCal Bioinformatics Inc. 2019
# 
# mass/charge helper functions

mass_proton <- function(){1.00727646688}
mass_isotope <- function(){1.0025}
mass_averagine <- function(){114.2446}

neutralMass <- function(mz, z=1){
  (mz - mass_proton()) * z
}

chargeMass <- function(mz, z=1){
  mz / z + mass_proton()
}
