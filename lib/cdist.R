# Jeff Jones
# SoCal Bioinformatics Inc. 2019
#
# cluster helper functions

# cartiesian distance
cdist <- Vectorize(function(a, b){
    sqrt(a^2+b^2)
})
