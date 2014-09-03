#if (!isGeneric('biodyn')) 
   setGeneric('biodyn', function(model,params,...)  standardGeneric('biodyn'))

utils::globalVariables('validParams')



# setParams<-function(model="pellat",its=1)
#   return(FLPar(as.numeric(NA),dimnames=list(params=c(validParams(model),"b0","q","sigma"),iter=its)))

   utils::globalVariables('validParams')
   models=factor(c('fox',      'schaefer',
                   'pellat',   'gulland',
                   'fletcher', 'shepherd',
                   'logistic', 'genfit'))
   
   modelParams=function(mdl) {
     if (is.factor(mdl)) mdl=as.character(mdl)
     list(fox       =c('r','k'),
          schaefer  =c('r','k'),
          pellat    =c('r','k','p'),
          shepherd  =c('r','k','m'),
          gulland   =c('r','k'),
          fletcher  =c('k','msy','p'),
          logistic  =c('k','msy'),
          genfit    =c('r','k','p'))[[mdl]]}
   
   defaultParams<-function(object) {
     params(object)<-FLPar(as.numeric(NA),dimnames=list(params=c(validParams(model(object)),'b0','q','sigma'),iter=1:dims(object)$iter))
     
     unt<-NULL
     if ('r'     %in% dimnames(params(object))$params){
       params(object)['r',    ]<-0.5
       unt<-c(unt,'')}
     if ('k'     %in% dimnames(params(object))$params){
       params(object)['k',    ]<-mean(catch(object))*10
       unt<-c(unt,units(catch(object)))}
     if ('p'     %in% dimnames(params(object))$params){
       params(object)['p',    ]<-2
       unt<-c(unt,'')}
     if ('msy'   %in% dimnames(params(object))$params){
       params(object)['msy',  ]<-mean(catch(object))
       unt<-c(unt,units(catch(object)))}
     if ('b0'    %in% dimnames(params(object))$params){
       params(object)['b0',   ]<-1
       unt<-c(unt,'')}
     if ('m'     %in% dimnames(params(object))$params){
       params(object)['m',    ]<-0.5
       unt<-c(unt,'')}
     if ('q'     %in% dimnames(params(object))$params){
       params(object)['q',    ]<-1.0
       unt<-c(unt,'')}
     if ('sigma' %in% dimnames(params(object))$params){
       params(object)['sigma',]<-0.3
       unt<-c(unt,'')}
     
     units(params(object))<-unt
     
     invisible(params(object))}
   
   # setParams<-function(model='pellat',its=1)
   #   return(FLPar(as.numeric(NA),dimnames=list(params=c(validParams(model),'b0','q','sigma'),iter=its)))
   
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

#" biodyn class
#"
#" @description A class to represent a biomass dynamic stock assessment model. 
#" There are slots for data, parameters and residuals and methods for calculating 
#" reference points and other derived quantities.
#" 
#" @return biodyn object
#" @export
#" @examples
#" \dontrun{biodyn()}
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
  params      =FLPar(c(.5,as.numeric(NA),2,1,as.numeric(NA),as.numeric(NA)),                            dimnames=list(params=c("r","k","p","b0"),iter=1)),
  control     =FLPar(array(rep(c(1,as.numeric(NA),as.numeric(NA),as.numeric(NA)),each=4), dim=c(4,4,1), dimnames=list(params=c("r","k","p","b0"),option=c("phase","min","val","max"),iter=1))),
  priors      =array(rep(c(0,0,0.3,1),       each=7), dim=c(7,4),   dimnames=list(params=c("r","k","p","b0","msy","bmsy","fmsy"),c("weight","a","b","type"))),
  vcov        =FLPar(array(as.numeric(NA), dim=c(4,4,1), dimnames=list(params=c("r","k","p","b0"),params=c("r","k","p","b0"),iter=1))),
  hessian     =FLPar(array(as.numeric(NA), dim=c(4,4,1), dimnames=list(params=c("r","k","p","b0"),params=c("r","k","p","b0"),iter=1))),
  objFn       =FLPar(array(as.numeric(NA),dim=c(2,1),dimnames=list("value"=c("ll","rss"),iter=1)))),
validity=validity) 

setGeneric("control",     function(object,method,...) standardGeneric("control"))
setGeneric("control<-",   function(object,value)      standardGeneric("control<-"))

setGeneric("diags",       function(object,method,...) standardGeneric("diags"))
setGeneric("diags<-",     function(object,value)      standardGeneric("diags<-"))

setMethod("catch",   signature(object="biodyn"),function(object, ...)   object@catch)
setMethod("catch<-", signature(object="biodyn", value="FLQuant"),
          function(object, value){
            updateFLComp(object, "catch", value)
            return(object)})

createFLAccesors <- function(class, exclude=character(1), include=missing) {
  
  object <- class

  if(!missing(include))
  	slots <- getSlots(class)[include]
  else
  	slots <- getSlots(class)[!names(getSlots(class))%in%exclude]

	defined <- list()

	for (x in names(slots)) {
		# check method is defined already and signatures match
		eval(
		substitute(if(isGeneric(x) && names(formals(x)) != "object") {warning(paste("Accesor
			method for", x, "conflicts with a differently defined generic. Type", x,
			"for more information")); break}, list(x=x))
			)
		# create new generic and accesor method
		eval(
		substitute(if(!isGeneric(x)) setGeneric(x, function(object, ...) standardGeneric(x)),
		list(x=x))
		)
		eval(
		substitute(setMethod(x, signature(y), function(object) return(slot(object, x))),
      list(x=x, y=class))
		)
		# create replacement method
		xr <- paste(x, "<-", sep="")
		eval(
		substitute(if(!isGeneric(x)) setGeneric(x,
			function(object, ..., value) standardGeneric(x)), list(x=xr))
		)
		eval(
		substitute(setMethod(x, signature(object=y, value=v), function(object, value)
			{slot(object, s) <- value; if(validObject(object)) object else stop("")}),
      list(x=xr, y=class, s=x, v=unname(slots[x])))
		)
    if(any(unname(slots[x]) %in% c("FLArray", "FLQuant", "FLCohort", "FLPar")))
    eval(
		substitute(setMethod(x, signature(object=y, value="numeric"), function(object, value)
			{slot(object, s)[] <- value; object}), list(x=xr, y=object, s=x))
		)
		defined[[x]] <- c(x, xr, paste("alias{",x,",", class,"-method}", sep=""),
			paste("\alias{",xr,",", class,",",unname(slots[x]), "-method}", sep=""),
			paste("\alias{",x,"-methods}", sep=""),
			paste("\alias{",xr, "-methods}", sep="")
		)
	}
	return(defined)
}	# }}}

invisible(createFLAccesors("biodyn", 
            exclude=c("desc","range","priors","objFn","mng","diags","stock",
                      "refpts","hessian",
                      "profile","mngVcov","ll"))) #,"priors","diags","objFn","control","mng")))
