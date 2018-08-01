library(fda)
library(fda.usc)
library(ncdf4)

setwd("Z://Flake")
load("newflakemat.Rdata", verbose=TRUE)
dim(newflakemat)


makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

makeTransparent2 <- function(someColor, alpha=100, levels=7)
{
  start <- 255-(levels*8)  
  newseq <- round(seq(start, 255, length.out=levels))
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=  newseq, maxColorValue=255)})
}

nc <- nc_open("Z://Flake/erai_tuned_runs//ALID399_oflake.nc")
time <- as.POSIXct(ncvar_get(nc,'time')*86400,origin='1970-01-01', tz ='UTC')
time1 <- strptime(time, format="%Y-%m-%d") #
dtime <- round((time1$year + time1$yday/365 + 1900),2)
foy <- time1$yday %/% 7 
foy

ex.nc <- nc_open("Z://Globolakes//ArcLakeClustering//ALID9999_CGREC9D_TS024LM.nc")
## to get time index
time2 <- as.POSIXct(ncvar_get(ex.nc,'time')*86400,origin='1970-01-01', tz ='UTC')
time22 <- strptime(time2, format="%Y-%m-%d") #
dtime2 <- round((time22$year + time22$yday/365 + 1900),2)
foy11 <- (time22$yday %/% 7) [1:398]
id <- ncvar_get(ex.nc, "lakeid")
idd <- which(id%in%rownames(newflakemat))
z <- (ncvar_get(ex.nc, "lswt")[idd,]-273.15)[,1:398]
dim(z)



all(id[idd]==rownames(newflakemat))

fulllam <- df2lambda(c(foy11), basis=create.bspline.basis(c(0,52), nbasis=53), df=12)


seasonalfd <- function(dat, wk){
  
  rounded <- round(unlist(lapply(split(dat, wk), mean),2))
  zeros= as.numeric(names(rounded[rounded<0.5]))+1
  shortBasis <- create.bspline.basis(c(0,52), nbasis=53)
  if (length(zeros)==0){
    ff <- Data2fd(dat, argvals=wk, basis=shortBasis, lambda=fulllam)
    nco <- ff$coefs}
  
  if (length(zeros)>0){
    mylam <- df2lambda(c(wk), basis=create.bspline.basis(c(0,52), nbasis=53-length(zeros)), df=12)
    redshortBasis <- create.bspline.basis(c(0,52), nbasis=53, dropind=zeros)
    ff <- Data2fd(dat, argvals=wk, basis=redshortBasis, lambda=mylam)
    nco <- rep(0, 53)
    nco[-zeros] <- ff$coefs
  }
  
  nco
}


origcoefs <- flakecoefs <- matrix(NA, nrow=732, ncol=53)
for (i in 1:732){
  bee <- seasonalfd(dat=z[i,], wk=foy11)
  cee <- seasonalfd(dat=newflakemat[i,], wk=foy)
  origcoefs[i,] <- bee
  flakecoefs[i,] <- cee
  print(i)
}

arcfd   <- fd( t(origcoefs), create.bspline.basis(c(0,52), nbasis=53))
flakefd <- fd( t(flakecoefs), create.bspline.basis(c(0,52), nbasis=53))

#save(arcfd, file="Z://Flake/May2017/arcfd.Rdata")
#save(flakefd, file="Z://Flake/May2017/flakefd.Rdata")
#basis <- create.bspline.basis(c(0,52), nbasis=53)
#save(basis, file="Z://Flake/May2017/basis.Rdata")


arc.fpc <- pca.fd(arcfd, 2)
cumsum(arc.fpc$varprop)
par(mfrow=c(1,2))
plot(arc.fpc)   
varmx.pca.fd(arc.fpc,2)
plot(arc.fpc$scores)

par(mfrow=c(1,1), cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
sc10 <- matrix(NA,ncol=2,nrow=732)
for(i in 1:732){
  sc10[i,] <- inprod(arc.fpc$harmonics, flakefd[i]-mean(flakefd))
  print(i)
}
plot(arc.fpc$scores, col=2, pch=20, xlab="FPC Score 1", ylab="FPC Score 2")
points(sc10, col=4, pch=20)
legend("topleft", col=c(2,4), pch=19, c("Arc", "FLake"), bty="n")


load("Z://CEHArcLake//11cols.Rdata", verbose=TRUE)
load("Z://Globolakes//ArcLakeClustering//globalcoefs398.Rdata", verbose=TRUE)


kkk=9
library(RColorBrewer)
set.seed(100); mycols <- sample(brewer.pal(12, "Set3"), 12)
sc1 <- arc.fpc$scores[,1]
sc2 <- arc.fpc$scores[,2]
kk <- kmeans(globalcoefs398, kkk, nstart=1000, iter.max=50)
newgroup <- kk$cluster
sort(table(newgroup))
dat <- data.frame(cbind(newgroup, sc1, sc2))
myqda = qda(newgroup ~ sc1 + sc2, data = dat)
class <- predict(myqda,method="looCV")
sum(class$class!=newgroup)
which(class$class!=newgroup)

id[idd][class$class!=newgroup]==rd$GloboLakes_ID[class$class!=newgroup]
class$class!=newgroup

write.csv(newrd, file="Z://Flake/May2017/lakeinfo_9groups_may2017.csv" )
