utils::globalVariables(c('calcR','fwdWindow','interval', 'nits',
                         'best', 'sim', 'i', 'rcvPeriod','pMeasure',
                         'dbWriteTable'))

#' runMSE
#' @description Runs a full MSE using an \code{FLStock} object as the Operating Model and \code{biodyn} as the Mangement Procedure
#'           
#' @aliases runMSE
#' 
#' @param om an \code{FLStock} objectl 
#' @param brp an \code{FLBRP} object that holds the biological parameters for use in the projections
#' @param range a \code{vector} the starting and end years for the projections, and the interval for running the MP
#' @param srDev  a \code{FLQuant} with recruitment deviates
#' @param uDev an \code{FLQuant} with stock recruitment residuals
#' @param hcr a \code{vector} with HCR options
#' @param ctrl  a \code{vector} with HCR options
#' @param con  a \code{vector} with HCR options
#' @param append  a \code{vector} with HCR options
#'
#'
#' @return invisibly a list of \code{data.frame}s with performance measures from OM and summaries from MP, if \code{con!=NULL} will
#' also write to a MYSQL database
#'  
#' @export
#' @rdname runMSE
#' 
#' @seealso \code{\link{biodyn}}
#' 
#' @examples
#' \dontrun{
#' simBiodyn()
#'    }
runMSE=function(om,brp,srDev,uDev,
                range=c(min=range(om)[max],max=range(om)[max]+30,interval=3),
                hcr  =c(ftar=0.75,fmin=0.01,blim=0.8,btrig=0.4),
                ctrl =controlFn(r=calcR(brp),k=refpts(brp)['msy','stock']*2),
                con=NULL,append=TRUE){
  
  ## OM projections
  om =fwdWindow(om,end=end+interval,brp)
  
  ## F in longterm
  lgt=FLQuant(refpts(brp)['msy','harvest']*hcr['ftar'],
                         dimnames=list(year=range['min']:(range['max']+range['interval']-1),iter=seq(nits)))
  
  ## Add stochastcity
 # om =fwd(om,f=fbar(FLCore::iter(om,1))[,ac(2:(range['min'])],sr=brp,sr.residuals=srDev)
  om =fwd(om,f=lgt, sr=brp,sr.residuals=srDev)
  
  ## save projection for comparison
  prj=om
  
  #### Set up MP
  ## SA Options
  ctrl=with(best[as.numeric(sim),3:6],controlFn(r,k,p,b0,phaseR=-1,phaseK=1))
  ctrl['r','phase']=hcr[i,'phaseR']
  ctrl['k','phase']=hcr[i,'phaseK']
  if (hcr[i,'p']!='known') ctrl['p',c('min','val','max')]=c(.1,1,10)
  
  res=mseBiodyn(om,brp,ctrl,srDev=srDev,
                         start=start+rcvPeriod,end=end,interval=3,
                         uCV =hcr[i,'uCV'],
                         ftar=hcr[i,'ftar'],blim=hcr[i,'blim'],btrig=hcr[i,'btrig'])
  
  prj=cbind(OM=sim,MP=i,pMeasure(res$prj,brp))
  mse=cbind(OM=sim,MP=i,pMeasure(res$om, brp))
  mp =cbind(OM=sim,MP=i,res$mp)
  
  dbWriteTable(con, 'prj', prj, append=append)
  dbWriteTable(con, 'mse', mse, append=append)
  dbWriteTable(con, 'mp',  mp,  append=append)
  
  invisible(res)}  
