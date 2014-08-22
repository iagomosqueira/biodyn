#' Creates a list of biodyn
#'
#' @description Creates a list of biodyn biomass dynamic model classes.
#' @return list of biodyn objects
#' @export
#' @examples
#' \dontrun{biodyns()}
setGeneric('biodyns', function(object, ...) standardGeneric('biodyns'))

setGeneric('biodyn',   function(model,params,...)  standardGeneric('biodyn'))

if (!isGeneric("msy"))     setGeneric('msy',      function(object,params,...) standardGeneric('msy'))
if (!isGeneric("fmsy"))    setGeneric('fmsy',     function(object,params,...) standardGeneric('fmsy'))
if (!isGeneric("bmsy"))    setGeneric('bmsy',     function(object,params,...) standardGeneric('bmsy'))
if (!isGeneric("refpts"))  setGeneric('refpts',   function(object,params,...) standardGeneric('refpts'))
if (!isGeneric("refptSE")) setGeneric('refptSE',  function(object,params,...) standardGeneric('refptSE'))

if (!isGeneric("harvest")) setGeneric('harvest',  function(object,params,...) standardGeneric('harvest'))

#if (!isGeneric("fwd"))     
  setGeneric("fwd",      function(object, ctrl, ...)    standardGeneric("fwd"))
#if (!isGeneric("hcr"))     
  setGeneric("hcr",      function(object, ...)          standardGeneric("hcr"))
#if (!isGeneric("hcrPlot")) 
  setGeneric("hcrPlot",  function(object, ...)          standardGeneric("hcrPlot"))
#if (!isGeneric("tac"))     
  setGeneric("tac",      function(object, harvest, ...) standardGeneric("tac"))

setGeneric('fit',   function(object,index,...)     standardGeneric('fit'))

if (!isGeneric("power"))    setGeneric('power',     function(object,ref,...)    standardGeneric('power'))
if (!isGeneric("diags"))    setGeneric('diags',     function(object,method,...) standardGeneric('diags'))
if (!isGeneric("diags<-"))  
  setGeneric('diags<-',   function(object,value)      standardGeneric('diags<-'))

if (!isGeneric("kobe")) 
   setGeneric('kobe',          function(object,method,...)    standardGeneric('kobe'))

#if (!isGeneric("control"))    
  setGeneric('control',     function(object,method,...) standardGeneric('control'))
#if (!isGeneric("control<-"))  
  setGeneric('control<-',   function(object,value)      standardGeneric('control<-'))

if (!isGeneric("plot"))  
  setGeneric('plot',   function(x,y)      standardGeneric('plot'))
  
#setGeneric('kobe',       function(object,method,...)    standardGeneric('kobe'))
if (!isGeneric("kobePhase")) setGeneric('kobePhase',  function(object,...)           standardGeneric('kobePhase'))
          
