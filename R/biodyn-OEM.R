#setMethod('survey', signature(object='FLStock'),
survey=function(object,timing=FLQuant(0,dimnames=dimnames(m(object))),wt=stock.wt(object),sd=0,log=FALSE,...){
        
      if(units(harvest(object)) == 'f'){
            res      <- colSums(stock.n(object)*exp(-harvest(object)*timing - m(object)*timing)*wt, na.rm=FALSE)
            dim(res) <- c(1, dim(res))
            dmns     <-dimnames(stock(object))
            dmns$iter<-dimnames(res)$iter
            
            res=FLQuant(res, dimnames=dmns)
      } else if(units(harvest(object)) == 'hr') {
            
            res      = colSums(object@stock.n*(1-object@harvest*timing)*exp(-object@m*timing)*timing*wt)
            dim(res) = c(1, dim(res))
            
            res=FLQuant(res, dimnames=dimnames(object@stock))
      } else
            stop("Correct units (f or hr) not specified in the harvest slot")
      
      if (sd>0){
           if (!log) res=apply(res,2:6,function(x,sd) rnorm( 1,x,sd=sd),sd=sd) else 
                     res=apply(res,2:6,function(x,sd) rlnorm(1,x,sd=sd),sd=sd)
           }
        
      res}

#setMethod('cpue', signature(object='FLStock'),
u=function(object,partialf=FLQuant(1,dimnames=dimnames(m(object))),wt=catch.wt(object),sd=0,log=FALSE,...){
            
            if(units(harvest(object)) == 'f'){
              res      <- colSums(stock.n(object)*harvest(object)*partialf/(harvest(object)+m(object))*exp(1-harvest(object)-m(object)), na.rm=FALSE)
              dim(res) <- c(1, dim(res))
              dmns     <-dimnames(stock(object))
              dmns$iter<-dimnames(res)$iter
              
              res=FLQuant(res, dimnames=dmns)
            } else if(units(harvest(object)) == 'hr') {
              
              res      = colSums(object@stock.n*(1-object@harvest*partialf)*exp(-object@m)*wt)
              dim(res) = c(1, dim(res))
              
              res=FLQuant(res, dimnames=dimnames(object@stock))
            } else
              stop("Correct units (f or hr) not specified in the harvest slot")
          
          if (sd>0){
            if (!log) res=apply(res,2:6,function(x,sd) rnorm( 1,x,sd=sd),sd=sd) else 
                      res=apply(res,2:6,function(x,sd) rlnorm(1,x,sd=sd),sd=sd)}
          
          res} 

cpueBiodym2Aspic=function(bd,type="B0",n=1,sd=0,log=TRUE){
  
   type=toupper(type)
   effort=harvest(bd)
   if (n>1)  effort=propagate(effort,n)
  
   dmns=dimnames(stock(bd))
   dmns$iter=1              
   dev=if (log) rlnorm(n,FLQuant(0,dimnames=dmns),sd) else rnorm( n,FLQuant(0,dimnames=dmns),sd)
  
   ctc=stock(bd)[,-dims(bd)$year]*effort
  
   switch(type,
          
          ## Error on Catch
          CC={ res=cbind(name="CC sim",model.frame(mcf(FLQuants(effort=effort,catch=ctc,dev=dev)),drop=T, stringsAsFactors=FALSE))
               res=transform(res,catch=catch*dev,index=catch*dev/effort)
               res},
          
          ## Error on index
          CE={ res=cbind(name="CE sim",model.frame(mcf(FLQuants(effort=effort,catch=ctc,dev=dev)),drop=T, stringsAsFactors=FALSE))
               res=transform(res,index=catch*dev/effort)
             },   
          
          ## Error on stock
          B0=cbind(name="B0 sim",model.frame(mcf(FLQuants(stock=stock(bd)*dev)),drop=T, stringsAsFactors=FALSE)),
          B1=cbind(name="I1 sim",model.frame(mcf(FLQuants(stock=(stock(bd)[,-dim(stock(bd))[2]]+stock(bd)[,-1]/2),dev=dev)),drop=T, stringsAsFactors=FALSE)),
          B2={ res=cbind(name="I0 sim",model.frame(mcf(FLQuants(stock=stock(bd),dev=dev)),drop=T, stringsAsFactors=FALSE))
               res},          
          
          ## Error on stock
          I0=cbind(name="B0 sim",model.frame(mcf(FLQuants(index=.1*stock(bd)*dev)),drop=T, stringsAsFactors=FALSE)),
          I1=cbind(name="I1 sim",model.frame(mcf(FLQuants(index=.1*(stock(bd)[,-dim(stock(bd))[2]]+stock(bd)[,-1]/2),dev=dev)),drop=T, stringsAsFactors=FALSE)),
          I2={ res=cbind(name="I0 sim",model.frame(mcf(FLQuants(index=.1*stock(bd),dev=dev)),drop=T, stringsAsFactors=FALSE))
               res}                
   )}



