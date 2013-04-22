combineFlr=function(...){

  object=list(...)
  if ("list" %in% is(object[[1]])) object=object[[1]]
  dmns=dimnames(object[[1]])
  dmns$iter=seq(sum(laply(object,function(x) dims(x)$iter)))
 
  if ("FLQuant" %in% is(object[[1]]))  res=FLQuant(NA,dimnames=dmns)
  if ("FLPar"   %in% is(object[[1]]))  res=FLPar(  NA,dimnames=dmns)
  
  its=llply(object,function(x) seq(dims(x)$iter))
  
  nits=0
  for (i in seq(length(its))){
     iter(res,its[[i]]+nits)=object[[i]]
     nits=max(its[[i]])}

  res}

combineFlr=function(...){
  
  object=list(...)
  if ("list" %in% is(object[[1]])) object=object[[1]]
  dmns=dimnames(object[[1]])
  dmns$iter=seq(sum(laply(object,function(x) dim(x)[length(dim(x))])))
  
  res=array(NA,unlist(laply(dmns,length)), dimnames=dmns)
  
  its=llply(object,function(x) dim(x)[length(dim(x))])
  
  nits=0
  for (i in seq(length(its))){
    res[,,,,,its[[i]]+nits]=object[[i]]
    nits=max(its[[i]])}
  
  res}