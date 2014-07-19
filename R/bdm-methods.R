# as.mcmc.rjags <- function (x) {
#   x <- x$BUGSoutput
#   if (x$n.chains > 1) {
#     z <- list()
#     for (i in 1:x$n.chains) {
#       z[[i]] <- mcmc(x$sims.array[, i, ], start = 1, thin = x$n.thin)
#     }
#     class(z) <- "mcmc.list"
#   }
#   else {
#     z <- mcmc(x$sims.matrix, start = 1, thin = x$n.thin)
#   }
#   return(z)
# }
# 
# setGeneric("k"                    ,function(object, ...)  standardGeneric("k"))            
# # setGeneric("r"                    ,function(object, ...)  standardGeneric("r"))            
# setGeneric("x"                    ,function(object, ...)  standardGeneric("x"))            
# setGeneric("q"                    ,function(object, ...)  standardGeneric("q"))            
# setGeneric("sigmaosq"             ,function(object, ...)  standardGeneric("sigmaosq"))            
# setGeneric("sigmapsq"             ,function(object, ...)  standardGeneric("sigmapsq"))            
# setGeneric("h"                    ,function(object, ...)  standardGeneric("h"))            
# # setGeneric("m"                    ,function(object, ...)  standardGeneric("m"))            
# setGeneric("g"                    ,function(object, ...)  standardGeneric("g"))             
# # setGeneric("stock"                ,function(object, ...)  standardGeneric("stock"))            
# setGeneric("depletion"            ,function(object, ...)  standardGeneric("depletion"))           
# # setGeneric("harvest"              ,function(object, ...)  standardGeneric("harvest_rate"))       
# setGeneric("sp"                   ,function(object, ...)  standardGeneric("sp")) 
# setGeneric("epsilonO"             ,function(object, ...)  standardGeneric("epsilon_o"))   
# setGeneric("epsilonP"             ,function(object, ...)  standardGeneric("epsilon_p"))  
# setGeneric("current_biomass"      ,function(object, ...)  standardGeneric("current_biomass"))    
# setGeneric("current_depletion"    ,function(object, ...)  standardGeneric("current_depletion"))  
# setGeneric("current_harvest_rate" ,function(object, ...)  standardGeneric("current_harvest_rate"))
# # setGeneric("bmsy"                 ,function(object, ...)  standardGeneric("bmsy"))
# # setGeneric("fmsy"                 ,function(object, ...)  standardGeneric("fmsy"))
# # setGeneric("msy"                  ,function(object, ...)  standardGeneric("msy"))
# # setGeneric("dmsy"                 ,function(object, ...)  standardGeneric("dmsy"))            
# # setGeneric("index"                ,function(object, ...)  standardGeneric("index"))
# # setGeneric("indexHat"             ,function(object, ...)  standardGeneric("indexHat"))
# # setGeneric("lp"                   ,function(object, ...)  standardGeneric("lp"))
#           
# setMethod("k"                    ,signature=c("bdm"), function(object) exp(object@trace$logK))             
# setMethod("r"                    ,signature=c("bdm","missing"), function(m) m@trace$r)             
# setMethod("x"                    ,signature=c("bdm"), function(object) object@trace$x)             
# setMethod("q"                    ,signature=c("bdm"), function(object) object@trace$q)            
# setMethod("sigmaosq"             ,signature=c("bdm"), function(object) object@trace$sigmaosq)             
# setMethod("sigmapsq"             ,signature=c("bdm"), function(object) object@trace$sigmaosq)            
# setMethod("h"                    ,signature=c("bdm"), function(object) object@trace$h)            
# setMethod("m"                    ,signature=c("bdm"), function(object) object@trace$m)             
# setMethod("g"                    ,signature=c("bdm"), function(object) object@trace$g)             
# setMethod("stock"                ,signature=c("bdm"), function(object) object@trace$biomass)             
# setMethod("depletion"            ,signature=c("bdm"), function(object) object@trace$depletion)           
# setMethod("harvest"              ,signature=c("bdm"), function(object) object@trace$harvest_rate)        
# setMethod("sp"                   ,signature=c("bdm"), function(object) object@trace$surplus_production)  
# setMethod("epsilonO"             ,signature=c("bdm"), function(object) object@trace$epsilon_o)    
# setMethod("epsilonP"             ,signature=c("bdm"), function(object) object@trace$epsilon_p)   
# setMethod("current_biomass"      ,signature=c("bdm"), function(object) object@trace$current_biomass)    
# setMethod("current_depletion"    ,signature=c("bdm"), function(object) object@trace$current_depletion)   
# setMethod("current_harvest_rate" ,signature=c("bdm"), function(object) object@trace$current_harvest_rate)
# setMethod("bmsy"                 ,signature=c("bdm"), function(object) object@trace$biomass_at_msy)
# setMethod("fmsy"                 ,signature=c("bdm"), function(object) object@trace$harvest_rate_at_msy)
# setMethod("msy"                  ,signature=c("bdm"), function(object) fmsy(object)/bmsy(object))
# setMethod("dmsy"                 ,signature=c("bdm"), function(object) object@trace$dmsy)            
# 
# setMethod("index"                ,signature=c("bdm"), function(object) object@trace$observed_index)
# setMethod("indexHat"             ,signature=c("bdm"), function(object) object@trace$predicted_index)
# setMethod("lp"                   ,signature=c("bdm"), function(object) object@trace$lp__)
