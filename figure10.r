source('default.R')

modelname=c("emilieTEST0427MR001minus17keepDD0518PLC3","11NOV0606gpsi1kf20copygmax1400day15")
biomassdata=rep(NA,8*length(modelname)*2)
dim(biomassdata)=c(8,length(modelname),2)
for (num in 2:length(modelname)) {
	model=modelname[num]
	for (ex in 1:2)  {
		new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"stomate.nc"))
		biomass=ncvar_get(new,'TOTAL_M')
		#mon=mymean(biomass,8)
		biomassdata[,num,ex]=toyear(biomass)/100
		
	}
}


obsBiomass=read.csv("caxBiomass.csv",header=F)

tiff("biomass 11June.tiff",width=600,height=400)
par(mar=c(5,5,3,3))
plot(2001:2008,biomassdata[,1,1],col='white',ylim=c(210,320),pch=19,xlab="",
ylab=expression(Biomass~density~(MgC~ha^-1)),cex.lab=2,cex.axis=2,cex=2,type='l',lwd=2)
ltys=c(1,2)
colos=c('black','red')
for (mo in 2:2) {
	for (ex in 1:2) {
		lines(2001:2008,biomassdata[,mo,ex],lwd=2,lty=3-ltys[mo],col=colos[ex])
	}
}


ctl_shift=obsBiomass[1:11,2]+(obsBiomass[12,2]-obsBiomass[1,2])
#points(obsBiomass[1:11,1],obsBiomass[1:11,2],pch=7,cex=2)
points(obsBiomass[1:11,1],ctl_shift,pch=7,cex=2)
points(obsBiomass[1:11,1],obsBiomass[12:22,2],pch=7,cex=2,col="red")
delta=obsBiomass[12,2]-obsBiomass[1,2]
lines(c(obsBiomass[1,1],obsBiomass[1,1]),
c(ctl_shift[1]+delta,ctl_shift[1]-delta),lwd=2)

legend("topright",c("CTL","TFE"),text.col=c("black","red"),bty="n",cex=2)
legend("topleft",legend=c("Model","Obs"),lty=c(1,-1),pch=c(-1,7),bty="n",cex=2)


#Arrows(x0 = 2008, y0 =biomassdata[8,2,1] , x1 = 2008, y1 = biomassdata[8,2,2],xpd=TRUE,lwd=2,arr.width=0.15,arr.type='simple')
#Arrows(x0 = 2008, y0 = biomassdata[8,2,2], x1 = 2008, y1 = biomassdata[8,2,1],xpd=TRUE,lwd=2,arr.width=0.15,arr.type='simple')


#Arrows(x0 = obsBiomass[7,1], y0 = ctl_shift[7], x1 = obsBiomass[7,1], y1 = obsBiomass[18,2],xpd=TRUE,lwd=2,arr.width=0.15,arr.type='simple')
#Arrows(x0 = obsBiomass[7,1], y0 = obsBiomass[18,2], x1 = obsBiomass[7,1], y1 = ctl_shift[7],xpd=TRUE,lwd=2,arr.width=0.15,arr.type='simple')

#text(2006.7,230,paste0(round(ctl_shift[7]-obsBiomass[18,2]),' MgC/ha'),cex=1.5)
#text(2007.2,270,paste0(round(biomassdata[8,2,1]-biomassdata[8,2,2]),' MgC/ha'),cex=1.5)

Arrows(x0 = 2002, y0 = 190, x1 = 2002, y1 = 205,xpd=TRUE,lwd=2,arr.width=0.15)

dev.off()

### 
tiff("biomass relative 9June.tiff",width=600,height=350)
par(mar=c(5,7,3,3))
plot(2001:2008,biomassdata[,1,1]/biomassdata[1,1,1]*100,col='white',ylim=c(83,103),pch=19,xlab="",
ylab="Biomass change \nrelative to 2001 (%)",cex.lab=2,cex.axis=2,cex=2,type='l',lwd=2,xaxt='n',yaxt='n')
ltys=c(1,2)
colos=c('black','red')
for (mo in 2:2) {
	for (ex in 1:2) {
		lines(2001:2008,biomassdata[,mo,ex]/biomassdata[1,mo,1]*100,lwd=2,lty=3-ltys[mo],col=colos[ex])
	}
}
ctl_shift=obsBiomass[1:11,2]+(obsBiomass[12,2]-obsBiomass[1,2])
#points(obsBiomass[1:11,1],obsBiomass[1:11,2],pch=7,cex=2)
points(obsBiomass[1:11,1],ctl_shift/ctl_shift[1]*100,pch=7,cex=2)
points(obsBiomass[1:11,1],obsBiomass[12:22,2]/ctl_shift[1]*100,pch=7,cex=2,col="red")
delta=obsBiomass[12,2]-obsBiomass[1,2]
#lines(c(obsBiomass[1,1],obsBiomass[1,1]),
#c(ctl_shift[1]+delta,ctl_shift[1]-delta),lwd=2)

legend("right",c("CTL","TFE"),text.col=c("black","red"),bty="n",cex=2)
legend("bottomleft",legend=c("Model","Obs"),lty=c(1,-1),pch=c(-1,7),bty="n",cex=2)

abline(h=100,lty=2,lwd=2)
axis(1,at=seq(2000,by=1,length.out=10),labels=seq(2000,by=1,length.out=10),cex.axis=2,lwd=2)
axis(2,at=seq(80,by=5,length.out=8),labels=seq(80,by=5,length.out=8),cex.axis=2,lwd=2)

axis(3,at=seq(2000,by=1,length.out=10),labels=FALSE,cex.axis=2,lwd=2,tck=0)
axis(4,at=seq(80,by=5,length.out=8),labels=FALSE,cex.axis=2,lwd=2,tck=0)

Arrows(x0 = 2002, y0 = 78, x1 = 2002, y1 = 81,xpd=TRUE,lwd=2,arr.width=0.15)

dev.off()
