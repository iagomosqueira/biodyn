utils::globalVariables('read.fit')

# io.ADMB.R - 
# FLCore/R/io.ADMB.R

# Copyright 2003-2012 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# $Id:  $

utils::globalVariables(c('admbCor','maply',
                         'str_trim',
                         'laply',
                         'llply',
                         'upper.triangle',
                         'str_trim'))


readADMB<-function(file){
  ## read in data from ADMB Par or Rep file
  dat  <-scan(file,what='',sep='\n',skip=1,quiet=TRUE)
  ## Get values
  vals <-lapply(strsplit(dat[grep('#',dat,invert=TRUE)],' '), function(x) as.numeric(x[nchar(x)>0]))
  ## name elements
  names(vals)<-lapply(grep('#',dat,value=T),function(x) substr(x,3,nchar(x)))

  return(vals)}

writeADMB<-function(x,file,append=FALSE){
  cat('#', names(x[1]),'\n',file=file,append=append)
  cat(x[[1]],'\n',file=file,append=TRUE)
  
  if (length(x)>1)
  for (i in 2:length(x)){
    cat('#', names(x[i]),'\n',file=file,append=TRUE)
    cat(x[[i]],'\n',file=file,append=TRUE)}}


read.admb <-
function(ifile)
{	
	fit			<- read.fit(ifile)
	class(fit)	<- c('list', 'admb')
	
	fn			<- paste(ifile,'.rep', sep='')
	rep			<- read.rep(fn)
	class(rep)	<- c('list', 'iscam')

	
	pfn			<- paste(ifile,'.psv',sep='')
	psv			<- NULL
	if(file.exists(pfn))
	{
		psv			<- read.psv(pfn)
		class(psv)	<- 'psv'
	}
	
	A			<- list()
	A$fit		<- fit
	A$rep		<- rep
	A$psv		<- psv
	return(A)
}

admbFit <-
function(ifile)
{
	# __Example:             
	#	file <-('~/admb/simple')
	#	A <- reptoRlist(file)
	#	Note there is no extension on the file name.
	
	## The following is a contribution from:
	## Anders Nielsen that reads the par & cor files.
	ret<-list() 
	parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
	 what='', n=16, quiet=TRUE)[c(6,11,16)]) 
	ret$nopar<-as.integer(parfile[1]) 
	ret$nlogl<-parfile[2] 
	ret$maxgrad<-parfile[3] 
	file<-paste(ifile,'.cor', sep='') 
	lin<-readLines(file) 
	ret$npar<-length(lin)-2 
	ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
	sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
	ret$names<-unlist(lapply(sublin,function(x)x[2])) 
	ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
	ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
	ret$cor<-matrix(NA, ret$npar, ret$npar) 
	corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
	ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
	ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
	ret$cov<-ret$cor*(ret$std%o%ret$std)
	
	class(ret) <-'admb'
	return(ret)
}

read.rep <- 
function(fn)
{
	# The following reads a report file
	# Then the 'A' object contains a list structure
	# with all the elemements in the report file.
	# In the REPORT_SECTION of the AMDB template use 
	# the following format to output objects:
	#  	report<<'object \n'<<object<<endl;
	#
	# The part in quotations becomes the list name.
	# Created By Steven Martell
	options(warn=-1)  #Suppress the NA message in the coercion to double
	
	
	ifile=scan(fn,what='character',flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
	idx=sapply(as.double(ifile),is.na)
	vnam=ifile[idx] #list names
	nv=length(vnam) #number of objects
	A=list()
	ir=0
	for(i in 1:nv)
	{
		ir=match(vnam[i],ifile)
		if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
		dum=NA
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=''))
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrows=irr-ir-1,fill=TRUE))

		if(is.numeric(dum))#Logical test to ensure dealing with numbers
		{
			A[[vnam[i]]]=dum
		}
	}
	options(warn=0)
	
	return(A)
}

read.psv <-function(fn, nsamples=10000){
	#This function reads the binary output from ADMB
	#-mcsave command line option.
	#fn = paste(ifile,'.psv',sep='')
	filen <- file(fn, 'rb')
	nopar <- readBin(filen, what = integer(), n = 1)
	mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
	mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
	close(filen)
	return(mcmc)}

getADMBHessian <- function(hess='admodel.hes'){
  ## This function reads in all of the information contained in the
  ## admodel.hes file. Some of this is needed for relaxing the covariance
  ## matrix, and others just need to be recorded and rewritten to file so ADMB
  ## 'sees' what it's expecting.
  filename <- file(hess, 'rb')
  on.exit(close(filename))
  num.pars <- readBin(filename, 'integer', 1)
  hes.vec <- readBin(filename, 'numeric', num.pars^2)
  hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
  hybrid_bounded_flag <- readBin(filename, 'integer', 1)
  scale <- readBin(filename, 'numeric', num.pars)
  result <- list(num.pars=num.pars, hes=hes,
                 hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
  return(result)
}

cv=function(hess){
  res         <- getADMBHessian(hess)
  cov         <- solve(res$hes[-1570,-1570])
  scale       <- res$scale
  cov.bounded <- cov*(scale %o% scale)
  se          <- sqrt(diag(cov.bounded))
  
  cov.bounded}

cv2=function(obj){
  cov         <- obj$hes
  scale       <- obj$scale
  cov.bounded <- cov*(scale %o% scale)
  se          <- sqrt(diag(cov.bounded))
  
  cov.bounded}

#file='/home/laurie/Desktop/gcode/mse4mfcl/ALB/papers/SCRS/SCRS2013-117/Inputs/pella.hst'

admbProfile=function(file) {
  
  dat =scan(file, what = '', sep = '\n', skip = 1,quiet=TRUE)
  nms =c('',unlist(lapply(grep('#', dat, value = T), function(x) str_trim(substr(x,2, nchar(x))))))
  
  brks=seq(length(dat))[substr(dat,1,1)=='#']
  brks=data.frame('from'=c(1,brks+1),'to'=c(1,brks[-1]-1,length(dat)))
  
  vals=mlply(brks, function(from,to,dat) as.numeric(unlist(strsplit(unlist(dat[from:to]),' '))),dat=dat)
  names(vals)=nms
  
  # when do profiles start?
  flag=seq(length(nms))[substr(nms,1,4)=='rand']
  
  prf=ldply(vals[-(1:(flag))], function(x) t(array(x,c(2,length(x)/2))))
  names(prf)=c('param','value','p')
  
  vals=vals[1:flag]
  vals[['profile']]=prf
  
  return(vals)}


admbCor<-function(ifile)
  {
    # __Example:             
    #	file <-('~/admb/simple')
    #	A <- reptoRlist(file)
    #	Note there is no extension on the file name.
    
    ## The following is a contribution from:
    ## Anders Nielsen that reads the par & cor files.
    ret<-list() 
    parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
                             what='', n=16, quiet=TRUE)[c(6,11,16)]) 
    ret$nopar<-as.integer(parfile[1]) 
    ret$nlogl<-parfile[2] 
    ret$maxgrad<-parfile[3] 
    file<-paste(ifile,'.cor', sep='') 
    lin<-readLines(file) 
    ret$npar<-length(lin)-2 
    ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
    sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
    ret$names<-unlist(lapply(sublin,function(x)x[2])) 
    ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
    ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
    ret$cor<-matrix(NA, ret$npar, ret$npar) 
    corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
    ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
    ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
    ret$cov<-ret$cor*(ret$std%o%ret$std)
    
    class(ret) <-'admb'
    return(ret)}


