source('default.R')

model="11NOV0606gpsi1kf20copygmax1400day15"
ind_show=rep(NA,8*2)
dim(ind_show)=c(8,2)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"stomate.nc"))
	ind=ncvar_get(new,"IND")
	ind_show[,ex]=toyear(ind)
}	


cck_show=rep(NA,8*2*20)
dim(cck_show)=c(8,2,20)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"stomate.nc"))
	for (co in 1:20) {
		ind=ncvar_get(new,paste0("CCK_",sprintf("%03d",co)))
		cck_show[,ex,co]=toyear(ind)*365
	}
}	


ccn_show=rep(NA,8*2*20)
dim(ccn_show)=c(8,2,20)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"stomate.nc"))
	for (co in 1:20) {
		ind=ncvar_get(new,paste0("CCN_",sprintf("%03d",co)))
		ccn_show[,ex,co]=ind[seq(1,2920,by=365)]
	}
}



ccd_show=rep(NA,8*2*20)
dim(ccd_show)=c(8,2,20)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"stomate.nc"))
	for (co in 1:20) {
		ind=ncvar_get(new,paste0("CCD_",sprintf("%03d",co)))
		ccd_show[,ex,co]=toyear(ind)
	}
}


weightmean=function(cck,ccn) {
	return(100*sum(cck)/sum(ccn))
}

bin_cck=rep(NA,8*2*4)
dim(bin_cck)=c(8,2,4)
for (year in 1:8) {
	for (ex in 1:2) {
		rows=which(ccd_show[year,ex,]<=0.2)
		bin_cck[year,ex,1]=weightmean(cck_show[year,ex,rows],ccn_show[year,ex,rows])
		rows=which(ccd_show[year,ex,]>0.2 & ccd_show[year,ex,]<=0.4)
		bin_cck[year,ex,2]=weightmean(cck_show[year,ex,rows],ccn_show[year,ex,rows])
		rows=which(ccd_show[year,ex,]>0.4)
		bin_cck[year,ex,3]=weightmean(cck_show[year,ex,rows],ccn_show[year,ex,rows])
		rows=which(ccd_show[year,ex,]>=0)
		bin_cck[year,ex,4]=weightmean(cck_show[year,ex,rows],ccn_show[year,ex,rows])
	}
}
load('/home/orchidee04/yyao/rebuild/daCosta.Rdata')

tiff("bin 003 11June.tiff",width=900,height=600)
part_name=c("<20cm (#1-#4)","20-40cm (#5-#6)",">40cm (#7-#20)","ALL (#1-#20)")
par(mfrow=c(2,2))
for (part in 1:4)  {
	par(mar=c(5,5,3,3))
	# open contrl gray TFE  ctl black 
	plot(bin_cck[,1,part],col="black",type="l",lty=1,lwd=2,
	main=part_name[part],cex.lab=2,cex.axis=2,
	cex.main=2,xaxt="n",xlab="",ylab=expression(Stem~mortality~rate~('%'~yr^-1)),
	ylim=c(0,10))
	#points(bin_cck[,1,part],pch=19,lwd=2)
	
	lines(bin_cck[,2,part],col="red",lwd=2)
	#points(bin_cck[,2,part],pch=19,col="red")
	axis(1, at=seq(0,by=1,length.out=11), 
					labels=seq(2000,by=1,length.out=11),cex.axis=2,lwd=2)
	axis(3, at=seq(0,by=1,length.out=11),labels=FALSE,cex.axis=2,lwd=2,tck=0)

	axis(2, at=seq(-2,12,by=2),labels=FALSE,lwd=2)
	axis(4, at=seq(-2,12,by=2),labels=FALSE,lwd=2,tck=0)
					
	if (part==1) {
		# <0.2m
		tfe_obs=mortality_b[1:8,1]
		ctl_obs=mortality_b[9:16,1]
		points(tfe_obs,pch=7,lwd=2,col="red")
		points(ctl_obs,pch=7,lwd=2)
		#lines(tfe_obs,lty=2,lwd=2,col="red")
		#lines(ctl_obs,lty=2,lwd=2)
	}
	if (part==2) {
		# <0.2m
		tfe_obs=mortality_c[1:8,1]
		ctl_obs=mortality_c[9:16,1]
		points(tfe_obs,pch=7,lwd=2,col="red")
		points(ctl_obs,pch=7,lwd=2)
		#lines(tfe_obs,lty=2,lwd=2,col="red")
		#lines(ctl_obs,lty=2,lwd=2)
	}
	if (part==3) {
		# <0.2m
		tfe_obs=mortality_d[1:8,1]
		ctl_obs=mortality_d[9:16,1]
		points(tfe_obs,pch=7,lwd=2,col="red")
		points(ctl_obs,pch=7,lwd=2)
		#lines(tfe_obs,lty=2,lwd=2,col="red")
		#lines(ctl_obs,lty=2,lwd=2)
	}
	if (part==4) {
		# <0.2m
		tfe_obs=mortality_a[1:8,1]
		ctl_obs=mortality_a[9:16,1]
		points(tfe_obs,pch=7,lwd=2,col="red")
		points(ctl_obs,pch=7,lwd=2)
		#lines(tfe_obs,lty=2,lwd=2,col="red")
		#lines(ctl_obs,lty=2,lwd=2)
	}
		
}

legend("topleft",c("CTL","TFE"),text.col=c("black","red"),bty="n",cex=2)
legend("topright",c("Model","Obs"),lty=c(1,-1),lwd=c(2,2),pch=c(-1,7),bty="n",cex=2)

dev.off()
