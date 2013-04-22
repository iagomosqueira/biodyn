
# biodyns {{{
biodyns <- setClass("biodyns", contains="FLComps",
                   validity=function(object) {
                     # All items are biodyn
                     if(!all(unlist(lapply(object, is, 'biodyn'))))
                       return("Components must be biodyn")	
                     
                     return(TRUE)})

setGeneric('biodyns', function(object, ...)
  standardGeneric('biodyns'))

# constructor
setMethod("biodyns", signature(object="biodyn"), function(object, ...) {
  lst <- c(object, list(...))
  biodyns(lst)
})

setMethod("biodyns", signature(object="missing"),
          function(...) {
            # empty
            if(missing(...)){
              new("biodyns")
              # or not
            } else {
              args <- list(...)
              object <- args[!names(args)%in%c('names', 'desc', 'lock')]
              args <- args[!names(args)%in%names(object)]
              do.call('biodyns',  c(list(object=object), args))
            }
          }
)

setMethod("biodyns", signature(object="list"),
          function(object, ...) {
            
            args <- list(...)
            
            # names in args, ... 
            if("names" %in% names(args)) {
              names <- args[['names']]
            } else {
              # ... or in object,
              if(!is.null(names(object))) {
                names <- names(object)
                # ... or in elements, ...
              } else {
                names <- unlist(lapply(object, name))
                # ... or 1:n
                idx <- names == "NA" | names == ""
                if(any(idx))
                  names[idx] <- as.character(length(names))[idx]
              }
            }
            
            # desc & lock
            args <- c(list(Class="biodyns", .Data=object, names=names),
                      args[!names(args)%in%'names'])
            
            return(
              do.call('new', args)
            )
            
          }) # }}}