setMethod('biodyn', signature(model="FLStock",params="missing"),
    function(model){
                  
      res      =new("biodyn")
      res@catch=catch(model)
      res@stock=window(stock(model),end=dims(model)$maxyear+1)
      
      dmns=dimnames(res@params)
      
      res@params=FLPar(as.numeric(NA),dimnames=dmns)
      res@params[]=c(.6,4*mean(res@catch),1,1)
      
      res@control[,"val"]=res@params
      res@control[,"min"]=res@control[,"val"]*.1
      res@control[,"max"]=res@control[,"val"]*10
      res@control[,"phase"]=c(1,1,-1,-1)
      
      range(res)[]=range(model)[c("minyear","maxyear")]
      
      res})


