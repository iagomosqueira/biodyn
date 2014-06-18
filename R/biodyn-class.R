utils::globalVariables("validParams")
models=factor(c("fox",      "schaefer",
                "pellat",   "gulland",
                "fletcher", "shepherd",
                "logistic", "genfit"))

modelParams=function(mdl) {
  if (is.factor(mdl)) mdl=as.character(mdl)
  list(fox       =c("r","k"),
       schaefer  =c("r","k"),
       pellat    =c("r","k","p"),
       shepherd  =c("r","k","m"),
       gulland   =c("r","k"),
       fletcher  =c("k","msy","p"),
       logistic  =c("k","msy"),
       genfit    =c("r","k","p"))[[mdl]]}

defaultParams<-function(object) {
  params(object)<-FLPar(NA,dimnames=list(params=c(validParams(model(object)),"b0","q","sigma"),iter=1:dims(object)$iter))
  
  unt<-NULL
  if ("r"     %in% dimnames(params(object))$params){
    params(object)["r",    ]<-0.5
    unt<-c(unt,"")}
  if ("k"     %in% dimnames(params(object))$params){
    params(object)["k",    ]<-mean(catch(object))*10
    unt<-c(unt,units(catch(object)))}
  if ("p"     %in% dimnames(params(object))$params){
    params(object)["p",    ]<-2
    unt<-c(unt,"")}
  if ("msy"   %in% dimnames(params(object))$params){
    params(object)["msy",  ]<-mean(catch(object))
    unt<-c(unt,units(catch(object)))}
  if ("b0"    %in% dimnames(params(object))$params){
    params(object)["b0",   ]<-1
    unt<-c(unt,"")}
  if ("m"     %in% dimnames(params(object))$params){
    params(object)["m",    ]<-0.5
    unt<-c(unt,"")}
  if ("q"     %in% dimnames(params(object))$params){
    params(object)["q",    ]<-1.0
    unt<-c(unt,"")}
  if ("sigma" %in% dimnames(params(object))$params){
    params(object)["sigma",]<-0.3
    unt<-c(unt,"")}
  
  units(params(object))<-unt
  
  invisible(params(object))}

# setParams<-function(model="pellat",its=1)
#   return(FLPar(NA,dimnames=list(params=c(validParams(model),"b0","q","sigma"),iter=its)))

getParams<-function(params,nm){
  if (nm %in% dimnames(params)$params)
    return(c(params[nm,]))
  else
    return(rep(as.numeric(NA),length=dims(params)$iter))}


validity<-function(object) {
  return(TRUE)
  ## Catch must be continous
  yrs<-dimnames(catch(object))$year
  
  if (!all(yrs == ac(dims(catch(object))$minyear:dims(catch(object))$maxyear)))
      return("years in catch not continous")

  # range
  dims <-dims(object)
  range<-as.list(object@range)

  return(TRUE)}

#' Biomass Dynamic Model Class
#'
#' @description A class that represents a biomass dynamic stock assessment model. It has slots
#' for data and parameters and methods for calculating reference points etc... Perversely there are 
#' no methods for estimating the parameters from data, the class needs to be extended to do this
#' , i.e. create a new class with estimation methods to estimate parameters
#' @return biodyn object
#' @export
#' @examples
#' \dontrun{aspic()}
 setClass("biodyn", representation(
    "FLComp",
    model         ="factor",
    catch         ="FLQuant",
    stock         ="FLQuant",
    diags         ="data.frame",
    params        ="FLPar",
    control       ="FLPar",
    priors        ="array",
    vcov          ="FLPar",
    hessian       ="FLPar",
    objFn         ="FLPar",
    ll            ="FLPar",
    mng           ="FLPar",
    mngVcov       ="FLPar",
    profile       ="data.frame"
    ),
  prototype(
    range       =unlist(list(minyear=as.numeric(NA), maxyear=as.numeric(NA))),
    catch       =FLQuant(),
    stock       =FLQuant(),
    model       =models[3],
    params      =FLPar(c(.5,NA,2,1,NA,NA),                            dimnames=list(params=c("r","k","p","b0"),iter=1)),
    control     =FLPar(array(rep(c(1,NA,NA,NA),each=4), dim=c(4,4,1), dimnames=list(params=c("r","k","p","b0"),option=c("phase","min","val","max"),iter=1))),
    priors      =array(rep(c(0,0,0.3,1),       each=7), dim=c(7,4),   dimnames=list(params=c("r","k","p","b0","msy","bmsy","fmsy"),c("weight","a","b","type"))),
    vcov        =FLPar(array(NA, dim=c(4,4,1), dimnames=list(params=c("r","k","p","b0"),params=c("r","k","p","b0"),iter=1))),
    hessian     =FLPar(array(NA, dim=c(4,4,1), dimnames=list(params=c("r","k","p","b0"),params=c("r","k","p","b0"),iter=1))),
    objFn       =FLPar(array(NA,dim=c(2,1),dimnames=list("value"=c("ll","rss"),iter=1)))),
	validity=validity) 