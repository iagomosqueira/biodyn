# setMethod("FLXSA", signature(stock="FLStock", indices="FLIndices"),
#           function(stock, indices, control=FLXSA.control(), diag.flag=TRUE) {
#             
#             Call <- match.call()
#             
#             # check FLIndices input
#             for (i in 1:length(indices)) {
#               # startf & endf present in @range
#               if (is.na(indices[[i]]@range["startf"]) || is.na(indices[[i]]@range["endf"]))
#                 stop(paste("Must supply startf & endf for range in FLIndex",i))
#               # @range has all elements
#               if (!all(c("min","max","plusgroup","minyear","maxyear","startf","endf") %in%
#                          names(indices[[i]]@range)))
#                 stop("Range must have names 'min','max','plusgroup','minyear','maxyear',
#                      'startf','endf'")
#               # adjust ranges to dims()
#               indices[[i]]@range["min"] <- max(indices[[i]]@range["min"], dims(indices[[i]])$min,
#                                                stock@range["min"])
#               indices[[i]]@range["max"] <- min(indices[[i]]@range["max"], dims(indices[[i]])$max,
#                                                stock@range["max"])
#               indices[[i]]@range["minyear"] <- max(indices[[i]]@range["minyear"],
#                                                    dims(indices[[i]])$minyear, stock@range["minyear"])
#               indices[[i]]@range["maxyear"] <- min(indices[[i]]@range["maxyear"],
#                                                    dims(indices[[i]])$maxyear, stock@range["maxyear"])
#               
#               # trim according to range
#               age <- indices[[i]]@range["min"]:indices[[i]]@range["max"]
#               year<- indices[[i]]@range["minyear"]:indices[[i]]@range["maxyear"]
#               
#               indices[[i]] <- trim(indices[[i]], age=age, year=year)
#             }
#             
#             # Double check validity
#             if (!validObject(stock))
#               stop("stock is not valid!")
#             if (!validObject(indices))
#               stop("FLIndices is not valid!")
#             if (!validObject(control))
#               stop("control is not valid!")
#             
#             # adjust range for minage
#             if ("minage" %in% names(stock@range))
#               minage <- stock@range["minage"]
#             else if ("min" %in% names(stock@range))
#               minage <- stock@range["min"]
#             else if ("minquant"  %in% names(stock@range))
#               minage <- stock@range["minquant"]
#             else
#               stop("'minage' not found in range")
#             
#             # adjust range for maxage
#             if ("maxage" %in% names(stock@range))
#               maxage <- stock@range["maxage"]
#             else if ("max" %in% names(stock@range))
#               maxage <- stock@range["max"]
#             else if ("maxquant" %in% names(stock@range))
#               maxage <- stock@range["maxquant"]
#             else
#               stop("'maxage' not found in range")
#             
#             # adjust plsugroup
#             if ("plusgroup" %in% names(stock@range))
#               stock@range["plusgroup"] <- maxage
#             
#             if (maxage<minage | stock@range["maxyear"] < stock@range["minyear"])
#               stop("Error in range")
#             
#             if (is.na(stock@range["plusgroup"]))
#               stop("Plus Group must be specified")
#             
#             # trim stock
#             stock <- trim(stock, year=stock@range["minyear"]:stock@range["maxyear"])
#             
#             if (!is.na(stock@range["plusgroup"]) &
#                   stock@range["plusgroup"] < dims(stock@catch.n)$max)
#               stock <- setPlusGroup(stock, stock@range["plusgroup"])
#             
#             stock@m <- stock@m[as.character(minage:maxage),,,,]
#             
#             # check catch.n is available
#             if (all(is.na(stock@catch.n)))
#               stop("catch.n is not available")
#             
#             stock@stock.n[] <- new('FLQuant')
#             stock@harvest[] <- new('FLQuant')
#             
#             # fqs
#             fqs <- function(assess) {
#               assess@index <- new("FLQuants", lapply(assess@index,FLQuant))
#               assess@index.hat <- new("FLQuants", lapply(assess@index.hat,FLQuant))
#               assess@index.var <- new("FLQuants", lapply(assess@index.var,FLQuant))
#               assess@index.res <- new("FLQuants", lapply(assess@index.res,FLQuant))
#               assess@q.hat <- new("FLQuants", lapply(assess@q.hat,FLQuant))
#               assess@q2.hat <- new("FLQuants", lapply(assess@q2.hat,FLQuant))
#               
#               if (validObject(assess))
#                 return(assess)
#               else
#                 stop("not valid")
#             }
#             
#             #
#             iters.stock  <-max(unlist(qapply(stock, function(x) dims(x)$iter)))
#             iters.indices<-max(unlist(lapply(indices@.Data,
#                                              function(x) max(unlist(qapply(x, function(x2) dims(x2)$iter))))))
#             
#             if ((iters.stock>1 && iters.indices>1) && missing(diag.flag))
#               diag.flag <- FALSE
#             
#             if ((iters.stock>1 && iters.indices>1) && diag.flag)
#               return("Multiple iters only allowed if diag.flag=FALSE")
#             
#             if(!diag.flag) {
#               
#               res<-.Call("FLXSA", iter(stock, 1), lapply(indices, iter, 1), control, FALSE)
#               iters <- max(iters.stock,iters.indices)
#               
#               if (iters>1) {
#                 
#                 res@stock.n<-propagate(FLQuant(res@stock.n@.Data),iters)
#                 res@harvest<-propagate(FLQuant(res@harvest@.Data),iters)
#                 for (i in as.character(2:iters)) {
#                   res. <- .Call("FLXSA", iter(stock,i), lapply(indices, iter,i), control, FALSE)
#                   iter(res@stock.n,i)<-FLQuant(res.@stock.n@.Data)
#                   iter(res@harvest,i)<-FLQuant(res.@harvest@.Data)
#                 }
#               }       
#               
#               res@harvest@units <- "f"
#               
#               return(res)
#             }
#             
#             res <-.Call("FLXSA", stock, indices, control, diag.flag)
#             
#             res <-fqs(.Call("FLXSA", stock, indices, control, diag.flag))
#             
#             if (class(res) != "FLXSA")
#               return(res)
#             res@call <- as.character(Call)
#             
#             # put wts amd nhats into a data frame
#             df <- as.data.frame(res@wts)
#             df1 <- (df[4])
#             df1[df1 >= 1, 1] <- paste("index", df1[df1 >= 1, 1])
#             df1[df1 == -1, 1] <- "fshk"
#             df1[df1 == -2, 1] <- "nshk"
#             df <- cbind(df[-4], df1)
#             
#             names(df) <- c("w", "nhat", "yrcls", "age", "year", "source")
#             
#             for(i in 1:length(indices)) {
#               v <- paste("index", i)
#               if (length(slot(indices[[i]],"name"))>0)
#                 df$source[df$source==v] <- slot(indices[[i]],"name")   
#             }
#             
#             wts.df <-df[,c(4,5,1,6)]
#             names(wts.df)[3] <-"data"
#             nhat.df<-df[,c(4,5,2,6)]
#             names(nhat.df)[3]<-"data"
#             
#             index.names<-sort(unique(df[,"source"]))
#             index.names<-index.names[substr(index.names,1,5)=="index"]
#             
#             fill.flq<-function(obj) {
#               dms <-dims(obj)
#               dmns<-dimnames(obj)
#               dmns[[1]]<-dms[[2]]:dms[[3]]
#               dmns[[2]]<-dms[[5]]:dms[[6]]
#               
#               res<-as.FLQuant(NA,dimnames=dmns)
#               res[dimnames(obj)[[1]],dimnames(obj)[[2]],,,]<-obj[dimnames(obj)[[1]],
#                                                                  dimnames(obj)[[2]],,,]
#               
#               return(res)
#             }
#             
#             res2  <- new("FLXSA")
#             res2@index.var<-new('FLQuants')
#             res2@index.hat<-new('FLQuants')
#             res2@index    <-new('FLQuants')
#             
#             j=0
#             
#             for (i in index.names) {
#               j=j+1
#               res2@index.var[[j]]<-1.0/fill.flq(as.FLQuant(wts.df[wts.df[,  "source"]==i,-4]))
#               res2@index.hat[[j]]<-    fill.flq(as.FLQuant(nhat.df[nhat.df[,"source"]==i,-4]))
#               res2@index.var[[j]]<-1.0/as.FLQuant(wts.df[wts.df[,  "source"]==i,-4])
#               res2@index.hat[[j]]<-as.FLQuant(nhat.df[nhat.df[,"source"]==i,-4])
#               dmns <- dimnames(res2@index.hat[[j]])
#               index  <- trim(indices[[j]]@index,age=dmns$age,year=dmns$year)
#               
#               res2@index[    is.na(index)]<-NA
#               res2@index.hat[is.na(index)]<-NA
#               res2@index.var[is.na(index)]<-NA
#               res2@index.res[is.na(index)]<-NA
#             }
#             
#             # F shrinkage
#             fshk  <- df[df[,"source"]=="fshk",]
#             
#             if (length(fshk[,1])>0) {
#               y.range <- range(fshk[,"year"])
#               a.range<-range(fshk[,"age"])
#               max.yr <-fshk[fshk[,"age"] ==a.range[2],]
#               max.age<-fshk[fshk[,"year"]==y.range[2],]
#               
#               res2@n.fshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
#                                                          year=y.range[1]:y.range[2]))
#               res2@n.fshk[as.character(max.age[,"age"]),as.character(y.range[2])]<-max.age[,
#                                                                                            "nhat"]
#               res2@n.fshk[as.character(a.range[2]),as.character(max.yr[,"year"])]<-max.yr[,"nhat"]
#               
#               res2@var.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
#                                                            year=y.range[1]:y.range[2]))
#               res2@var.nshk[as.character(max.age[,"age"]),as.character(y.range[2])]<-
#                 1/max.age[,"w"]
#               res2@var.nshk[ as.character(a.range[2]),as.character(max.yr[,"year"])]<-
#                 1/max.yr[,"w"]
#             }
#             
#             # N shrinkage
#             nshk   <-df[df[,"source"]=="nshk",]
#             if (length(nshk[,1])>0) {
#               y.range<-range(nshk[,"year"])
#               a.range<-range(nshk[,"age"])
#               max.yr <-nshk[nshk[,"age"] ==a.range[2],]
#               max.age<-nshk[nshk[,"year"]==y.range[2],]
#               
#               res2@n.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
#                                                          year=y.range[1]:y.range[2]))
#               res2@n.nshk[,as.character(y.range[2])] <- max.age[,"nhat"]
#               res2@n.nshk[ as.character(a.range[2])] <- max.yr[,"nhat"]
#               
#               res2@var.nshk <- as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
#                                                            year=y.range[1]:y.range[2]))
#               res2@var.nshk[,as.character(y.range[2])]<-1/max.age[,"w"]
#               res2@var.nshk[ as.character(a.range[2])]<-1/max.yr[,"w"]
#             }
#             
#             res2@diagnostics<- df
#             res2@index.hat  <- res@index.hat
#             res2@stock.n    <- FLQuant(res@stock.n@.Data)
#             res2@harvest    <- FLQuant(res@harvest@.Data)
#             res2@survivors  <- FLQuant(res@survivors@.Data)
#             res2@se.int     <- FLQuant(res@se.int@.Data)
#             res2@se.ext     <- FLQuant(res@se.ext@.Data)
#             res2@n.fshk     <- FLQuant(res@n.fshk@.Data)
#             res2@n.nshk     <- FLQuant(res@n.nshk@.Data)
#             res2@var.fshk   <- FLQuant(res@var.fshk@.Data)
#             res2@var.nshk   <- FLQuant(res@var.nshk@.Data)
#             res2@harvest@units <- "f"
#             res2@q.hat      <- res@q.hat
#             res2@q2.hat     <- res@q2.hat
#             res2@index.res  <- res@index.res
#             res2@index      <- res@index
#             
#             res2@index.hat  <-FLQuants(res@index.hat)
#             res2@index.res  <-FLQuants(res@index.res)
#             res2@index.var  <-FLQuants(res@index.var)
#             res2@q.hat      <-FLQuants(res@q.hat)    
#             res2@q2.hat     <-FLQuants(res@q2.hat)   
#             
#             for (i in 1:length(indices)) {
#               res2@index.range[[i]]<-indices[[i]]@range
#               res2@index.name[i]   <-indices[[i]]@name
#               res2@index[[i]]      <-res@index[[i]]
#               res2@index.hat[[i]]  <-FLQuant(res@index.hat[[i]],
#                                              dimnames=dimnames(res@index.hat[[i]]))
#               res2@index.res[[i]]  <-FLQuant(res@index.res[[i]],
#                                              dimnames=dimnames(res@index.res[[i]]))
#               res2@index.var[[i]]  <-FLQuant(res@index.var[[i]],
#                                              dimnames=dimnames(res@index.var[[i]]))
#               res2@q.hat[[i]]      <-FLQuant(res@q.hat[[i]],
#                                              dimnames=dimnames(res@q.hat[[i]]))
#               res2@q2.hat[[i]]     <-FLQuant(res@q2.hat[[i]],
#                                              dimnames=dimnames(res@q2.hat[[i]]))
#               
#               dmns<-dimnames(indices[[i]]@index)
#               na. <-is.na(indices[[i]]@index[
#                 dmns$age[dmns$age %in% dimnames(res2@index[[i]])$age],
#                 dmns$year[dmns$year %in% dimnames(res2@index[[i]])$year]])
#               
#               res2@index[[i]][na.] <- NA
#               res2@index.hat[[i]][na.]<- NA
#               res2@index.res[[i]][na.]<- NA
#               res2@index.var[[i]][na.]<- NA
#               
#               dimnames(res2@index[[i]]    )$unit   <-dimnames(indices[[i]]@index)$unit
#               dimnames(res2@index.hat[[i]])$unit   <-dimnames(indices[[i]]@index)$unit
#               dimnames(res2@index.res[[i]])$unit   <-dimnames(indices[[i]]@index)$unit
#               dimnames(res2@index.var[[i]])$unit   <-dimnames(indices[[i]]@index)$unit
#               dimnames(res2@q.hat[[i]])$unit       <-dimnames(indices[[i]]@index)$unit
#               dimnames(res2@q2.hat[[i]])$unit      <-dimnames(indices[[i]]@index)$unit
#               dimnames(res2@index[[i]]    )$season <-dimnames(indices[[i]]@index)$season
#               dimnames(res2@index.hat[[i]])$season <-dimnames(indices[[i]]@index)$season
#               dimnames(res2@index.res[[i]])$season <-dimnames(indices[[i]]@index)$season
#               dimnames(res2@index.var[[i]])$season <-dimnames(indices[[i]]@index)$season
#               dimnames(res2@q.hat[[i]])$season     <-dimnames(indices[[i]]@index)$season
#               dimnames(res2@q2.hat[[i]])$season    <-dimnames(indices[[i]]@index)$season
#               dimnames(res2@index[[i]]    )$area   <-dimnames(indices[[i]]@index)$area
#               dimnames(res2@index.hat[[i]])$area   <-dimnames(indices[[i]]@index)$area
#               dimnames(res2@index.res[[i]])$area   <-dimnames(indices[[i]]@index)$area
#               dimnames(res2@index.var[[i]])$area   <-dimnames(indices[[i]]@index)$area
#               dimnames(res2@q.hat[[i]])$area       <-dimnames(indices[[i]]@index)$area
#               dimnames(res2@q2.hat[[i]])$area      <-dimnames(indices[[i]]@index)$area   
#             }
#             
#             res2@control <- res@control
#             res2@call    <- res@call
#             res2@desc    <- res@desc
#             res2@range   <- stock@range
#             units(res2@harvest)<-"f"
#             
#             res2@control@shk.n   <-as.logical(res2@control@shk.n)
#             res2@control@shk.f   <-as.logical(res2@control@shk.f)
#             res2@control@vpa     <-as.logical(res2@control@vpa)
#             
#             return(res2)
#           }
#           )
# # }}}
#         