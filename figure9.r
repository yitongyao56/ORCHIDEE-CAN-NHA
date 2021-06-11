source('default.R')
mortality_annual=rep(NA,2*2*8*20)
dim(mortality_annual)=c(2,2,8,20) 

modelname="11NOV0606gpsi1kf20copygmax1400day15"
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname,"stomate.nc"))
	for (co in 1:20) {
		mortality=ncvar_get(new,paste0("circ_class_mor_",sprintf("%03d",co)))
		for (year in 1:8) {
			aa=mortality[((year-1)*365+1):(year*365)]
			mortality_annual[ex,1,year,co]=sum(aa)
		}
		ccn=ncvar_get(new,paste0("CCN_",sprintf("%03d",co)))
		for (year in 1:8) {
			aa=ccn[((year-1)*365+1)]
			mortality_annual[ex,2,year,co]=mortality_annual[ex,1,year,co]/aa
		}
	}
}

circ_class_mortality=rep(NA,2920*2*20)
dim(circ_class_mortality)=c(2920,2,20)

for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname,"stomate.nc"))
	for (co in 1:20) {
		mortality=ncvar_get(new,paste0("circ_class_mor_",sprintf("%03d",co)))
		circ_class_mortality[,ex,co]=mortality
	}
}



tiff("figure 9 psi 11June.tiff",width=1200,height=700)
#var_label=c("leaf storage","stem storage","root storage","psi_leaf","psi_stem","psi_root","PLC (%)")
par(mfrow=c(2,2))
for (ex in 1:2) {
		num=4
		par(mar=c(5,5,3,3))
		df=vl_daily_data[,ex,num,]
		colos=brewer.pal(8, "Spectral")
		if (num ==4) {
			colos=colos[8:1]
		}
		df[df>90]=89.9
		df[df<10]=10.1
		image.plot(list(x=1:2920,y=seq(0,by=1/20,length.out=20),z=df[,20:1]),las=1,zlim=c(10,90),
		col=colos,xaxt="n",ylab="Cohorts",cex.axis=2,breaks=c(10,20,30,40,50,60,70,80,90),
		cex.main=3,yaxt="n",legend.cex=2,cex.lab=2,legend.width=2,main=expe[ex],
		axis.args=list(cex.axis=2),yaxs="i",legend.shrink=1,legend.mar=12)
		axis(1,at=seq(-365,by=365,length.out=10),labels=seq(2000,by=1,length.out=10),cex.axis=2,lwd=2)
		axis(2,at=seq(-1/10,by=1/10,length.out=11),labels=c(' ',seq(20,by=-2,length.out=10)),cex.axis=2,lwd=2)
		axis(3,at=seq(-365,by=365,length.out=10),labels=FALSE,lwd=2,tck=0)
		axis(4,at=seq(-1/10,by=1/10,length.out=11),labels=FALSE,lwd=2,tck=0)
		mtext("PLC (%)",side=3,cex=2,line=1,adj=1.12)

	
	par(mar=c(5,5,3,3))
	colos=brewer.pal(7, "YlOrRd")
	#colos=c(colos)
	
	df=mortality_annual[ex,2,,]
	df[df==0]=-999
	#df[df>100]=0
	#df[df>70]=70
	df[df>0.16]=0.159
	image.plot(list(x=2001:2009,y=seq(0,by=1/20,length.out=20),z=df[,20:1]),
	las=1,zlim=c(0,0.16),col=colos[1:6],
	ylab="Cohorts",cex.axis=2,breaks=c(0,0.01,0.03,0.06,0.09,0.12,0.16),
	lab.breaks=c(0,0.01,0.03,0.06,0.09,0.12,0.16),xaxt='n',
	cex.main=3,legend.cex=2,cex.lab=2,legend.width=2,yaxt="n",main=expe[ex],yaxs="i",
	axis.args=list(cex.axis=2),legend.shrink=1,legend.mar=12)
	axis(1,at=seq(2000,by=1,length.out=10),labels=seq(2000,by=1,length.out=10),cex.axis=2,lwd=2)
	axis(2,at=seq(-1/10,by=1/10,length.out=12),labels=c(' ',seq(20,by=-2,length.out=11)),cex.axis=2,lwd=2)
	axis(3,,at=seq(2000,by=1,length.out=10),labels=FALSE,lwd=2,tck=0)
	axis(4,at=seq(-1/10,by=1/10,length.out=12),labels=FALSE,lwd=2,tck=0)
	mtext("Fraction",side=3,cex=2,line=1,adj=1.12)
}

dev.off()
