#' tseries
#' @description A utility method to create a data.frame with time series of Performance Measures i.e.
#' stock biomass, SSB, recruitment, F, harvest rate
#'            
#' @param stk    \code{object} 
#' @param refpts \code{FLBRP}
#' @param proxy  \code{character} with name of reference point in refpts(brp), by default 'msy'.
#'
#' @return a \code{data.frame} 
#'  
#' @aliases tseries,FLStock,missing  
#'  
#' @docType methods
#' @rdname tseries
#'  
#' @examples
#'  \dontrun{tseries(ple4,FLBRP(ple4))}
#'    
setGeneric('tseries',function(object,refpts,...)    standardGeneric('tseries'))

tseriesFn1=function(stk){
  
  res=FLQuants(stock     =stock(stk)%/%bmsy(stk),
               catch     =catch(stk)%/%msy( stk),
               harvest   =harvest(stk)%/%fmsy(stk))
  
  model.frame(res,drop=T)}

tseriesFn2=function(stk,brp,proxy='msy'){
  
  res=FLQuants(stock  =stock(stk)%/%refpts(brp)[proxy,'biomass'],
               ssb    =ssb(  stk)%/%refpts(brp)[proxy,'ssb'],
               rec    =rec(  stk)%/%refpts(brp)[proxy,'rec'],
               catch  =catch(stk)%/%refpts(brp)[proxy,'yield'],
               fbar   =fbar( stk)%/%refpts(brp)[proxy,'harvest'],
               harvest=(catch(stk)/stock(stk))%/%(refpts(brp)[proxy,'yield']/refpts(brp)[proxy,'biomass']))
  
  model.frame(res,drop=T)}
# 
# setMethod('tseries', signature(object='FLStock',refpts='FLBRP'), 
#   function(object,refpts,proxy='msy',...) {
# 
#   res=tseriesFn2(object,refpts,proxy)
#     
#   return(res)})

setMethod('tseries', signature(object='FLStock',refpts='missing'), 
          function(object,refpts,...) {
            
            res=tseriesFn1(object)
            
            return(res)})
