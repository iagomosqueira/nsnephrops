################################################################################
#                                                                              #
# Albacore N Atlantic MSE                                                      #
# MSE                                                                          #
#                                                                              #
################################################################################

#### Initial Stuff #############################################################
## FLR Packages
library(FLCore)
library(FLash)
library(FLBRP)
library(FLXSA)

## Observation Error Model
OEM<-function(object,start,type="stock",end=NULL,plusgroup=NULL,...)
    {
    if (missing(end)) end<-start
    yrs<-as.character(end:start)
    stk<-window(object, start=start, end=end)

    ## replace any slots
    args <- list(...)
    slt<-names(getSlots("FLStock"))[getSlots("FLStock")=="FLQuant"]
    for(i in names(args)[names(args) %in% names(slt)[slt=="FLQuant"]]){
       yrs.      <-yrs[yrs %in% dimnames(slot(res, i))$year]
       slot(res, i)[,yrs.]<-i[,yrs.]}

    if (!is.null(plusgroup))
       stk<-setPlusGroup(stk,plusgroup)

    return(stk)}

## CPUE Index
OEMIndex<-function(stk,start,end=NULL,plusgroup=NULL,startf=NULL,endf=NULL,deviates=NULL){
     ## unbiased population estimates
     if (missing(end)) end<-start
     yrs<-start:end

     stk<-window(stk,start=start,end=end)
     if (!is.null(plusgroup))
        stk<-setPlusGroup(stk,plusgroup)@n

     idx<-as(stk,"FLIndex")

     if (!is.null(startf)) idx@range["startf"]<-startf
     if (!is.null(endf))   idx@range["endf"]  <-endf

     if (!is.null(deviates))
        idx@index<-idx@index*deviates[dimnames(idx@index)$age,ac(yrs)]

     return(idx)}

### cas simulator
casSetup<-function(flq,Linf=314.9,K=0.0892,t0=-1.132,CV=0.3,niters=100){
   la       <-function(age,Linf,K,t0) Linf*(1-exp(-K*(age-t0)))

   dmns     <-dimnames(flq)
   dmns$unit<-1:niters
   cas      <-FLQuant(la(as.numeric(dmns$age),Linf,K,t0),dimnames=dmns)
   cas      <-cas*rlnorm(prod(dim(cas)),0,CV)

   return(cas)}
   
## function to apply a linearly increasing trend to an FLQuant
biasLinear<-function(x,obj){
   if (x>0)
      res  <-1-(sort(cumsum(rep(x, dims(obj)$year)),d=T)-x)
   else
      res  <-sort(1-(sort(cumsum(rep(x, dims(obj)$year)),d=T)-x),d=T)

   return(obj*FLQuant(rep(res,each=dims(obj)$age),dimnames=dimnames(obj)))}

hcrF<-function(SSB,Bpa,Blim,Fmin,Fmax){
    val <-qmax(qmin((Fmax-(Fmax-Fmin)*(Bpa-SSB)/(Bpa-Blim)),Fmax),Fmin)

    return(val)}

runMSE<-function(OM,srPar,srRsdl=srRsdl,fmult=1){
  nits<-dims(OM)$iter
   
  dmns  <-list(age=1:15,year=1970:2030,iter=1:nits)
  idxDev<-FLQuant(rlnorm( prod(unlist(lapply(dmns,length))),0,.3),dimnames=dmns)

  nits<-dims(OM)$iter
  iYr <-2005
  XSActrl<-FLXSA.control(maxit=80,tsrange=100,tspower=0)

  #### OEM setup
  MPstk<-OEM(     OM,start=1970,end=2029)
  MPidx<-OEMIndex(OM,start=1970,end=2029,startf=0,endf=0.25,deviates=idxDev)

  for (iYr in 2005:2028){
     cat("===================", iYr, "===================\n")
     #### OEM
     MPstk          <-window(MPstk,end  =iYr)
     MPidx          <-window(MPidx,end  =iYr)
     MPstk[,ac(iYr)]<-OEM(     OM, start=iYr)
     MPidx[,ac(iYr)]<-OEMIndex(OM, start=iYr,deviates=idxDev)

     #### Assessment
     MPstk      <-MPstk+FLXSA(MPstk,MPidx,XSActrl,diag.flag=FALSE)
     if (iYr==2005)
        MPbrp      <-brp(FLBRP(MPstk))

     #### Calculate TAC
     MP  <-stf(MPstk,nyears=2)
     ctrl<-fwdControl(data.frame(year=1:2+iYr,quantity=c("f","f")))
     dms<-dimnames(ctrl@trgtArray)
     dms$iter<-1:nits
     ctrl@trgtArray<-array(NA,lapply(dms,length),dms)
     ctrl@trgtArray[1,"val", ]<-mean(fbar(MP)[,ac(iYr-(1:3)),drop=T])
     ctrl@trgtArray[2,"val", ]<-MPbrp@refpts["msy","harvest",,drop=T]*fmult

     SRrs<-FLQuant(c(apply(rec(MP)[,ac(1970:(iYr-1))],6,function(x) exp(mean(log(x))))),dimnames=list(year=0:2+iYr,iter=1:nits))
     MP  <-fwd(MP,ctrl=ctrl,sr=list(model="mean",params=FLPar(1)),sr.residuals=SRrs)

     #### Go Fish
     TAC<-catch(MP)[,ac(iYr+2),drop=T]
     ctrl<-fwdControl(data.frame(year=iYr+2,max=c(NA,2),quantity=c("catch","f")))
     dms<-dimnames(ctrl@trgtArray)
     dms$iter<-1:nits
     ctrl@trgtArray<-array(NA,lapply(dms,length),dms)
     ctrl@trgtArray[1,"val", ]<-TAC
     ctrl@trgtArray[2,"max", ]<-2.0
     OM <-fwd(OM,ctrl=ctrl,sr=list(model="bevholt",params=srPar),sr.residuals=srRsdl)

     print(plot(lapply(FLStocks(OM,MP),window,end=iYr+2)))
     }

   plot(plot(lapply(FLStocks(OM,MP),window,end=iYr+2)))

   invisible(list(OM=OM,MP=MP))
   }

runMSEPoor<-function(OM,srPar,srRsdl,MPcas,refYrs=1975:1985,catchMult=1.0,idxMult=0.25){
  nits<-dims(OM)$iter
  XSActrl<-FLXSA.control(maxit=40,tsrange=100,tspower=0)

  #### OEM setup
  MPstk<-OEM(     OM,   start=1970,end=2029)
  MPcas<-window(  MPcas,start=1970,end=2029)
  MPIdx<-apply(sweep(MPcas,c(1:2,6),stock.n(OM)[,dimnames(MPcas)$year],"*"),c(2,6),mean)/10e6

  IdxTarget  <-apply(MPIdx[,ac(refYrs)],6,mean)
  IdxLimit   <-IdxTarget*idxMult
  catchTarget<-apply(catch(MPstk)[,ac(refYrs)],6,mean)*catchMult
  catchMax   <-catchTarget
  
  for (iYr in 2005:2029){
     cat("===================", iYr, "===================\n")

     #### OEM
     MPstk<-OEM(OM,iYr)
     MPIdx[,ac(iYr)]<-apply(sweep(MPcas[,ac(iYr)],c(1:2,6),stock.n(OM)[,ac(iYr)],"*"),c(2,6),mean)/10e6

     #### Calculate TAC
     IdxMean<-apply(MPIdx[,ac(-(5:1)+iYr)],6,mean)
     TAC<-qmin(catchTarget*qmax((IdxMean-IdxLimit)/(IdxTarget-IdxLimit),0.1), catchMax)

     #### Go Fish
     ctrl<-fwdControl(data.frame(year=iYr+1,max=c(NA,2),quantity=c("catch","f")))
     dms<-dimnames(ctrl@trgtArray)
     dms$iter<-1:nits
     ctrl@trgtArray<-array(NA,lapply(dms,length),dms)
     ctrl@trgtArray[1,"val", ]<-TAC
     ctrl@trgtArray[2,"max", ]<-2.0
     OM <-fwd(OM,ctrl=ctrl,sr=list(model="bevholt",params=srPar),sr.residuals=srRsdl)

     print(plot(lapply(FLStocks(OM),window,end=iYr+2)))
     }

   print(plot(lapply(FLStocks(OM),window,end=iYr+2)))

   invisible(list(OM=OM,MP=MPIdx))
   }


#load("\\\\Laurielpt/Stuff/FLR/flr4mse/data/albMSE.RData")
##load("C:/Stuff/FLR/flr4mse/data/albMSE.RData")
load("albMSE.RData")


dmns     <-dimnames(srRsdl)
dmns$iter<-1:100
dmns$year<-1970:2030
rho      <-0.3
srRsdl   <-FLQuant(sample(c(srRsdl),prod(unlist(lapply(dmns,length))),TRUE),dimnames=dmns)
for (i in as.numeric(dmns$year[-1]))
   srRsdl[,ac(i)]<-srRsdl[,ac(i)]*(1.0-rho)+srRsdl[,ac(i-1)]*rho

OM=stf(propagate(window(alb4B,start=1970),100),nyear=22)
harvest(OM)[,ac(2005:2030)]<-iter(harvest(OM),1)[,sample(ac(1996:2005),prod(dim(harvest(OM)[1,ac(2005:2030)])),TRUE)]

fmsy         <-list()
fmsy[["0.75"]]<-runMSE(OM, srPar4B, srRsdl=srRsdl,fmult=0.75)
fmsy[["0.5" ]]<-runMSE(OM, srPar4B, srRsdl=srRsdl,fmult=0.5)
fmsy[["1.0" ]]<-runMSE(OM, srPar4B, srRsdl=srRsdl,fmult=1.0)

resIdx         <-list()
resIdx[["1.00"]]<-list()
resIdx[["0.75"]]<-list()

MPcas=casSetup(m(OM))
resIdx[["0.75"]][["1"]]<-runMSEPoor(OM,srPar4B,srRsdl,MPcas,catchMult=0.75)
resIdx[["1.00"]][["1"]]<-runMSEPoor(OM,srPar4B,srRsdl,MPcas,catchMult=1.00)
resIdx[["1.25"]][["1"]]<-runMSEPoor(OM,srPar4B,srRsdl,MPcas,catchMult=1.25)

MPcas=casSetup(m(OM),CV=0.6,niters=50)
resIdx[["0.75"]][["2"]]<-runMSEPoor(OM,srPar4B,srRsdl,MPcas,catchMult=0.75)
resIdx[["1.00"]][["2"]]<-runMSEPoor(OM,srPar4B,srRsdl,MPcas,catchMult=1.00)
resIdx[["1.25"]][["2"]]<-runMSEPoor(OM,srPar4B,srRsdl,MPcas,catchMult=1.25)

plot(FLStocks(resIdx[["1.25"]][["1"]][["OM"]],resIdx[["1.25"]][["2"]][["OM"]],
              resIdx[["1.00"]][["1"]][["OM"]],resIdx[["1.00"]][["2"]][["OM"]],
              resIdx[["0.75"]][["1"]][["OM"]],resIdx[["0.75"]][["2"]][["OM"]],
              resVPA[["OM"]],resMSY[["OM"]]))
              
tmp<-FLStocks(resIdx[["1.25"]][["1"]][["OM"]],
              resIdx[["1.00"]][["1"]][["OM"]],
              resIdx[["0.75"]][["1"]][["OM"]])
save(tmp,file="c:/temp/ind1.RData")

tmp<-FLStocks(resIdx[["1.25"]][["2"]][["OM"]],
              resIdx[["1.00"]][["2"]][["OM"]],
              resIdx[["0.75"]][["2"]][["OM"]])
save(tmp,file="c:/temp/ind2.RData")

save(FLStocks(resIdx[["1.25"]][["1"]][["OM"]],resIdx[["1.25"]][["2"]][["OM"]],
              resIdx[["1.00"]][["1"]][["OM"]],resIdx[["1.00"]][["2"]][["OM"]],
              resIdx[["0.75"]][["1"]][["OM"]],resIdx[["0.75"]][["2"]][["OM"]],
              resVPA[["OM"]],resMSY[["OM"]]),file="c:/temp/MSE.RData")



brp<-FLBRP(fmsy[["0.75"]][["OM"]],sr=list(model="bevholt",param=srPar4B)))
brp<-FLBRP(resMSY[["OM"]],sr=srPar4B)