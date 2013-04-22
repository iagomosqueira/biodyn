#### ADMB ###################################################################################
## 1) setExe copies exe from bin and creates a temo dir
## 2) Write data files
## 3) Run exe
## 4) Read output files
#############################################################################################

#tmp=model.frame(FLPar(array(NA,dim=c(3,2,1),dimnames=list(param=c("a","b","c"),var=c("x","y"),iter=1))))
#tmp=as.data.frame(FLPar(array(NA,dim=c(3,2,1),dimnames=list(param=c("a","b","c"),var=c("x","y"),iter=1))))

## copies exe into temp dir ready to run
setExe=function(exeNm,package,dir=tempdir()){
  ##### set up temp dir with exe for data files
  # Linux
  if (R.version$os=="linux-gnu") {
    exe = paste(system.file("bin", "linux", package=package, mustWork=TRUE),exeNm, sep="/")
    file.copy(exe, dir)
    dir = paste(dir, "/", sep="")
    
    # Windows
  } else if (.Platform$OS.type == "windows") {
    exe = paste(system.file("bin", "windows", package=package, mustWork=TRUE), paste(exeNm, ".exe", sep=""), sep="/")
    file.copy(exe, dir)
    dir = paste(dir, "\\", sep="")
    
    # Mac OSX
  }else 
    stop()
  
  oldwd = getwd()
  
  # change wd to avoid exe case bug
  setwd(dir)
  
  oldwd}

setPella=function(obj, exeNm="pella", dir=tempdir()) {
  # create input files ################################
  dgts=options()$digits
  options(digits=22)

  # cpue
  if (is.FLQuant(obj[[2]])){  
     idx=model.frame(FLQuants(index=obj[[2]]), drop=TRUE)
     idx=data.frame(idx,name=1)
  }else if ("FLQuants" %in% class(obj[[2]])){  
     idx=as.data.frame(obj[[2]], drop=TRUE)
     names(idx)[2:3]=c("index","name")
     idx=transform(idx,name=as.numeric(name))
  }
  idx=idx[!is.na(idx$index),]

  bd.        =obj[[1]]
 
  nms=c( biodyn:::modelParams("pellat"),"b0")
  if (length(unique(idx$name))>0)
    nmIdx=paste(c("q","sigma"), rep(unique(idx$name),each=length(unique(idx$name))),sep="")
  else  
    nmIdx=c("q","sigma")

  
  # ctl file
  ctl        = bd.@control[nms,]
  ctl[,2:4]  = ctl[,c(2,4,3)]
  ctl        = alply(ctl,1)
  names(ctl) = nms

  biodyn:::writeADMB(ctl, paste(dir, "/", exeNm, ".ctl", sep=""),FALSE)

  cat("# q ####################\n", file=paste(dir, "/", exeNm, ".ctl", sep=""),append=TRUE)

  ctl           = bd.@control[nmIdx[grep("q",nmIdx)],]
  ctl[,2:4]     = ctl[,c(2,4,3)]
  ctl           = alply(t(matrix(ctl,dim(ctl))),1)
  names(ctl)    = c("phase","lower","upper","guess")

  biodyn:::writeADMB(ctl, paste(dir, "/", exeNm, ".ctl", sep=""),TRUE)
  
  cat("# sigma ################\n", file=paste(dir, "/", exeNm, ".ctl", sep=""),append=TRUE)
  ctl           = bd.@control[nmIdx[grep("s",nmIdx)],]
  ctl[,2:4]     = ctl[,c(2,4,3)]
  ctl           = alply(t(matrix(ctl,dim(ctl))),1)
  names(ctl)    = c("phase","lower","upper","guess")
  biodyn:::writeADMB(ctl, paste(dir, "/", exeNm, ".ctl", sep=""),TRUE)
   
  # prr file
  prr = bd.@priors[c(nms,nmIdx),] 
  prr = alply(prr,1)
  names(prr) = dimnames(bd.@priors)$params
  writeADMB(prr, paste(dir, "/", exeNm, ".prr", sep=""))
   
  # write data
  ctc = as.list(model.frame(bd.[["catch"]], drop=TRUE))
  ctc = c(nYrs=length(ctc[[1]]), ctc)
  res = c(ctc, c(nIdxYrs=dim(idx)[1], nIdx=length(unique(idx$name)), idx))
  
  writeADMB(res, paste(dir, "/", exeNm, ".dat", sep=""))
  
#   # propagate as required
#   its = dims(bd)$iter
#   
#   # params
#   bd@params = propagate(params(bd)[,1], its)
#   # stock
#   stock(bd)  = FLQuant(dimnames=dimnames(stock(bd))[1:5], iter=its)
#   
  # vcov
  vcov(bd.)=FLPar(array(NA, dim     =c(dim(params(bd.))[1],dim(params(bd.))[1],1), 
                           dimnames=list(params=dimnames(params(bd.))[[1]],
                                         params=dimnames(params(bd.))[[1]],iter=1)))
  
  
  options(digits=dgts)
  return(bd.)}

getPella=function(obj, exeNm="pella") {
  t1 = read.table(paste(exeNm,".rep",sep=""),skip =18,header=T) 

  # params
  t2 = unlist(c(read.table(paste(exeNm,".rep",sep=""),nrows=4)))
  q. = unlist(c(read.table(paste(exeNm,".rep",sep=""),nrows=1,skip=8)))
  s. = unlist(c(read.table(paste(exeNm,".rep",sep=""),nrows=1,skip=10)))

  nms=c("r","k","b0","p")
  obj@params[nms,] = t2
  
  obj@params[grep("q",dimnames(obj@params)$params),]=q. 
  obj@params[grep("s",dimnames(obj@params)$params),]=s. 
  t3 = unlist(c(read.table(paste(exeNm,".rep",sep=""),skip=dim(params(obj))[1]*2,nrows=2,header=F)))

  obj@objFn["ll"] =t3[length(t3)]
  obj@objFn["rss"]=t3[length(t3)-1]
  
  us=paste("u",seq(length(dimnames(params(obj))$params[grep("q",dimnames(params(obj))$params)])),sep="")
  obj@ll=FLPar(biodyn:::readADMB("lls.txt"),dimnames=list(params=us,iter=1))
  
  # stock biomass
  obj@stock[,1:dim(t1)[1]] = unlist(c(t1["stock"])) 
  
  return(obj)} 

#FLParBug sim@control["r","val",1]=c(.5,.6)
#exe(object, exeNm="biodyn", dir=tempdir(), set=biodyn:::set, get=biodyn:::set, cmdOps=paste("-maxfn 500"))
  
activeParams=function(obj) dimnames(control(obj))$params[c(control(obj)[,"phase"]>-1)]

## runs exe

##############################################################
#' fit
#'
#' Estimates parameters in a \code{biodyn} class by fitting catch to CPUE indices
#' 
#'
#' @param   \code{object}, an object of class \code{biodyn}
#' @param   \code{index}, an \code{FLQuant}, \code{FLQuants} or  \code{data.frame} object with CPUE indices
#' @param   \code{cmdOps}, a character string giving ADMB options see \url{http://www.admb-project.org/documentation/manuals/ADMBrefcard-A4.pdf/view}
#'
#' @export
#' @docType methods
#' @rdname fit
#'
#' @examples
#' /dontrun{
#' data(bd)
#' bd=fit(bd,swonIndex)}
setMethod("fit",signature(object='biodyn',index="FLQuant"),
          function(object,index=index,exeNm="pella",package="biodyn", 
                   dir=tempdir(),
                   set=setPella,
                   get=getPella,cmdOps=paste("-maxfn 500 -iprint 0"))
   fitPella(object,index=index,exeNm="pella",package="biodyn", 
            dir=tempdir(),
            set=setPella,
            get=getPella,cmdOps=paste("-maxfn 500 -iprint 0")))

setMethod("fit",signature(object='biodyn',index="FLQuants"),
          function(object,index=index,exeNm="pella",package="biodyn", 
                   dir=tempdir(),
                   set=setPella,
                   get=getPella,cmdOps=paste("-maxfn 500 -iprint 0"))
            fitPella(object,index=index,exeNm="pella",package="biodyn", 
                     dir=tempdir(),
                     set=setPella,
                     get=getPella,cmdOps=paste("-maxfn 500 -iprint 0")))

fitPella=function(object,index=index,exeNm="pella",package="biodyn", 
                  dir=tempdir(),
                  set=setPella,
                  get=getPella,cmdOps=paste("-maxfn 500 -iprint 0"))          
  {
  first=TRUE          
  if (dims(object)$iter==1 &  1<ifelse(is(index)[1]=="FLQuant",dims(index)$iter>1,max(laply(index,function(x) dims(x)$iter))))
    catch(object)=propagate(catch(object),dims(index)$iter)
  
  max=min(dims(catch(object))$maxyear,ifelse(is(index)[1]=="FLQuant",dims(index)$maxyear,max(laply(index,function(x) dims(x)$maxyear))))
  if (!is.na(range(object)["maxyear"])) max=min(max,range(object)["maxyear"])
  min=min(dims(catch(object))$minyear,ifelse(is(index)[1]=="FLQuant",dims(index)$minyear,max(laply(index,function(x) dims(x)$minyear))))
  if (!is.na(range(object)["minyear"])) min=max(min,range(object)["minyear"])

  object=window(object,start=min,end=max)
  if ("FLQuant" %in% is(index)){ 
    index =window(index,start=min,end=max)
  }else if ("FLQuants" %in% is(index)){
    index =FLQuants(llply(index, window,start=min,end=max))}
  
  slts=getSlots("biodyn")
  slts=slts[slts %in% c("FLPar","FLQuant")]

  #oldwd =setExe(exeNm,package,dir)
  oldwd=getwd()
  setwd(dir)
  exe()

  object=list(object,index)
  bd =object[[1]]
  its=max(maply(names(slts), function(x) dims(slot(bd,x))$iter))
  nms=dimnames(params(bd))$params
  bd@vcov   =FLPar(array(NA, dim=c(length(nms),length(nms),its), dimnames=list(params=nms,params=nms,iter=seq(its))))
  bd@hessian=bd@vcov
  
  us=paste("u",seq(length(dimnames(params(bd))$params[grep("q",dimnames(params(bd))$params)])),sep="")
  bd@ll=FLPar(NA,dimnames=list(params=us,iter=seq(1)))

  if (its>1){
      ## these are all results, so doesnt loose anything
      bd@stock  =iter(bd@stock,  1)
      bd@params =iter(bd@params, 1)
      bd@objFn  =iter(bd@objFn,  1)
      bd@vcov   =iter(bd@vcov,   1)
      bd@ll     =iter(bd@ll,     1)
      bd@hessian=iter(bd@hessian,1)
      bd@mng    =FLPar(a=1)
      
      bd=propagate(bd,its)
      }
   
  cpue=object[[2]]
  for (i in seq(its)){     
     object[[2]] = iter(cpue,i) 
    
     for (s in names(slts)[-(7:8)]){      
        slot(object[[1]],s) = iter(slot(bd,s),i) 
        }  
     
     object[[1]]=set(object,exeNm,dir)

     # run
     #system(paste("./", exeNm, " ", cmdOps, sep=""))
     system(paste(exeNm, " ", cmdOps, sep=""))
     
     # gets results
     object[[1]]=getPella(object[[1]], exeNm)
     
     for (s in names(slts)[slts=="FLQuant"]){
         iter(slot(bd,s),i) = slot(object[[1]],s)
         } 
     
     if (its<=1){
       ##hessian
       x<-file(paste(dir,"admodel.hes",sep="/"),'rb')
       nopar<-readBin(x,"integer",1)
       H<-matrix(readBin(x,"numeric",nopar*nopar),nopar)
       try(bd@hessian@.Data[activeParams(object[[1]]),activeParams(object[[1]]),i] <- H, silent=TRUE)
       close(x)
     
       ## vcov
       if (file.exists(paste(dir,"admodel.cov",sep="/")))
         try(bd@vcov@.Data[activeParams(object[[1]]),activeParams(object[[1]]),i] <- cv(), silent=TRUE) 
       #if (file.exists(paste(dir,"admodel.cov",sep="/"))){
       #   x<-file(paste(dir,"admodel.cov",sep="/"),'rb')
       #   nopar<-readBin(x,"integer",1)
       #   H<-matrix(readBin(x,"numeric",nopar*nopar),nopar)
       #   try(bd@vcov@.Data[activeParams(object[[1]]),activeParams(object[[1]]),i] <- H, silent=TRUE)
       #close(x)}
       }
   
     bd@params@.Data[  ,i] = object[[1]]@params
     bd@control@.Data[,,i] = object[[1]]@control
     #bd@objFn@.Data[   ,i] = object[[1]]@objFn

     #FLParBug
     bd@ll@.Data[,i][] = unlist(c(object[[1]]@ll))
     
     if (file.exists("pella.std")){
     
     err=try(mng.<-read.table("pella.std",header=T)[,-1])
             
     ## FLPar hack
     if (first & any(is(err)!="try-error")) {
       bd@mng   =FLPar(array(unlist(c(mng.[,-1])),dim=c(dim(mng.)[1],2,its),dimnames=list(param=mng.[,1],var=c("hat","sd"),iter=seq(its))))
       first=!first
    }else{
       if (is(err)!="try-error") bd@mng@.Data[,,i][]=unlist(c(mng.[,-1]))
       }}
  }
  
  units(bd@mng)="NA"
  
  bd=fwd(bd,catch=catch(bd)[,rev(dimnames(catch(bd))$year)[1]])
  
  if (length(grep("-mcmc",cmdOps))>0 & length(grep("-mcsave",cmdOps))>0){
    #"-mcmc 100000 -mcsave 100"
    setMCMC=function(obj,dir){
      ps=exp(read.psv(paste(dir,"pella.psv",sep="/")))
      ps=FLPar(t(ps))
      dimnames(ps)$params=activeParams(bd)
      units(ps)=NA
      ps}
    
    par=setMCMC(bd,dir)
    bd@params=propagate(bd@params[,1],dims(par)$iter)
    bd@params[dims(par)$params,]=par
    bd@stock=propagate(bd@stock,dim(params(bd))[2])
    bd=fwd(bd,catch=catch(t.))
    }
  
  if (its<=1) bd@diags=getDiags()
  
  setwd(oldwd)
                                  
  return(bd)}
            
#library(matrixcalc)
#is.positive.definite(H)

#vc=vcov(bd)[activeParams(bd),activeParams(bd),drop=T]
#vc[lower.triangle(vc)==0]=t(vc)[lower.triangle(vc)==0]
#is.positive.definite(vc,tol=1e-8)
#(vc-t(vc))/vc


admbCor=function(fl="pella.cor"){
  if (!file.exists(fl)) return
  require(matrixcalc)

  tmp=options()
  options(warn=-1)
  
  res=scan(fl,sep="\n",what=as.character())[-1]
  res=maply(res, str_trim)
  res=strsplit(res, " +")
  
  nms=laply(res,function(x) x[2])[-1]
  hat=laply(res,function(x) as.numeric(x[3]))[-1]
  sd =laply(res,function(x) as.numeric(x[4]))[-1]
  
  npar=length(res)-1
  
  val =as.numeric(unlist(llply(res[-1], function(x) strsplit(x," +")[-(1:4)])))
  cor=upper.triangle(array(1,dim=c(npar,npar),dimnames=list(nms,nms)))
  cor[upper.triangle(cor)==1]=val
  cor=cor+t(cor)
  diag(cor)=1
  
  vcov=cor*sd%o%sd
  
  options(tmp)
  
  hat =FLPar(hat, units="NA")
  dimnames(hat)$params=dimnames(vcov)[[1]]
  vcov=FLPar(vcov,units="NA")
  names(vcov)[1:2]=c("params","params")
  
  
  return(FLPars(hat =hat,
                vcov=vcov))}

getDiags=function(fl="pella.rep"){
  skip=grep("Model",scan(fl,sep="\n",what=as.character()))
  
  res=read.table(fl,skip=skip,header=TRUE)
  res=transform(res,residual=log(index/hat))
  res=cpue:::diagsFn(res)
  
  res$harvest=res$catch/res$stock
  res$stock  =res$stock.
  
  res[,-7]}

calcElasticity=function(bd,mn=3,rg=5){
  require(numDeriv)
  
  elasFn=function(x,dmns,bd,mn,rg) {
    
    params(bd)[dmns]=exp(x)
    bd=fwd(bd,catch=catch(bd))
    
    maxYr =ac(range(bd)["maxyear"])
    max1  =ac(range(bd)["maxyear"]-(seq(mn)-1))
    max2  =ac(range(bd)["maxyear"]-(seq(mn)-1+mn))
    maxR  =ac(range(bd)["maxyear"]-(seq(rg)-1))
     
    smy=c(  stock(bd)[,maxYr],
          harvest(bd)[,maxYr],
              msy(bd),
             fmsy(bd),
             bmsy(bd),
            catch(bd)[,maxYr]/msy( bd),
            stock(bd)[,maxYr]/bmsy(bd),
            harvest(bd)[,maxYr]/fmsy(bd),
            mean(stock(bd)[,max1])/mean(stock(bd)[,max2]),
            mean(harvest(bd)[,max1])/mean(harvest(bd)[,max2]),
            abs(coefficients(lm(data~year,as.data.frame(stock(  bd)[,maxR],drop=T)))["year"]),
            abs(coefficients(lm(data~year,as.data.frame(harvest(bd)[,maxR],drop=T)))["year"]))     
      
    return(log(smy))}
  
  parNms=c(biodyn:::modelParams(model(bd)),"b0")

  jbn=jacobian(elasFn,log(c(bd@params[parNms])),dmns=parNms,bd=bd,mn=mn,rg=rg)
  
  
  dimnms=list(params=parNms,
              stat  =c("stock",    "harvest",
                       "msy",      "bmsy",     "fmsy",
                       "catchMsy", "stockBmsy","harvestFmsy",
                       "stockMn",  "harvestMn",
                       "stockRg",  "harvestRg"))
  
  jbn=FLPar(array(t(jbn),dim=dim(t(jbn)),dimnames=dimnms))
  units(jbn)="NA"
  
  return(jbn)}

exe=function(package="biodyn"){

  sep =  function() if (R.version$os=="linux-gnu") ":" else if (.Platform$OS=="windows") ";" else ","
  
  # Linux
  if (R.version$os=="linux-gnu") {
    path = paste(system.file("bin", "linux",  package=package, mustWork=TRUE), sep="/")
    # Windows
  } else if (.Platform$OS.type == "windows") {
    path = paste(system.file("bin", "windows", package=package, mustWork=TRUE), sep="/")
    # Mac OSX
  }else 
    stop()
  
  #path="C:\\R\\R-2.15.1\\library\\aspic\\bin\\windows"
  
  path <- paste(path, sep(), Sys.getenv("PATH"),sep="")
   
  Sys.setenv(PATH=path)

return(path)}
#asp=aspic("http://gbyp-sam.googlecode.com//svn//trunk//tests//aspic//swon//2009//high//aspic.inp")
#asp=fit(asp)

calcSS=function(x) daply(x@diags, .(name), 
                         with, sum(residual^2,na.rm=T)/sum(count(!is.na(residual))))

#FLParBug
#tst=FLPar(NA,dimnames=list(params=c("a","b"),col=c("x","y"),iter=1:2))
#tst["a","x",2]=3


