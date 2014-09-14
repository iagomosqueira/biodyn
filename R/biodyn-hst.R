# fl  ="/tmp/RtmpzT3zmG/pella.hst"
# hst =scan(fl,sep="\n",what=as.character())
# hash=substr(hst,1,1)=="#"
# 
# pos  =seq(length(t.))[hash]
# start=pos+1
# end  =c(pos[-1]-1,length(hst))
# 
# names=hst[seq(length(t.))[hash]]
# names=str_trim(substr(names,2,nchar(names)))
# 
# key=data.frame(name=names,start=start,end=end)
# 
# array(unlist(strsplit(hst[22:53]," ")),c(2,32))

admbPlt<-function(fl){
  hst=scan(fl,sep="\n",what=as.character())
  if (length(hst)<10) return(NULL)  
  
  lns=(seq(length(hst))[substr(hst,1,7)=="Profile"]+1):
      (seq(length(hst))[substr(hst,1,7)=="Minimum"][1]-1)
  
  prfl=as.data.frame(t(array(as.numeric(unlist(strsplit(str_trim(hst[lns])," "))),
                             c(2,length(lns)))))
  
  lns=(seq(length(hst))[substr(hst,1,6)=="Normal"]+1):
      (seq(length(hst))[substr(hst,1,7)=="Minimum"][2]-1)
  
  nrm =as.data.frame(t(array(as.numeric(unlist(strsplit(str_trim(hst[lns])," "))),
                             c(2,length(lns)))))
  
  res =data.frame(prfl,nrm[,2])
  names(res)=c("value","p","pNorm")
  
  return(res)}
#fl ="/tmp/RtmpzT3zmG/lpr.plt"   
# prfl=mdply(data.frame(var=c("r", "k","bnow","fnow",
#                             "msy","bmsy","fmsy","cmsy",
#                             "bbmsy","ffmsy","bk","fr",
#                             "bratio","fratio","slopeb","slopef")),
#            function(var){
#                fl=paste("/tmp/RtmpzT3zmG/lp",var,".plt",sep="")
#                admbPlt(fl)})
# 
# prfl=ddply(prfl,.(var), transform, sp=cumsum(p)/sum(p))
# 
# ggplot(prfl)+geom_path(aes(value,sp))+facet_wrap(~var,scales="free")
