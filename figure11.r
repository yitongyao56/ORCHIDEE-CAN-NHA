source('default.R')

expe=c("CTL","TFE")
modelname="11NOV0606gpsi1kf20copygmax1400day15"
expe_data=rep(NA,2920*2*7)
dim(expe_data)=c(2920,2,7)

expe_sum=rep(NA,8*2*7)
dim(expe_sum)=c(8,2,7)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname,"stomate.nc"))   # 19maint
	leaf=ncvar_get(new,"BM_ALLOC_LEAF")
	sap_ab=ncvar_get(new,"BM_ALLOC_SAP_AB")
	sap_be=ncvar_get(new,"BM_ALLOC_SAP_BE")
	fruit=ncvar_get(new,"BM_ALLOC_FRUIT")
	root=ncvar_get(new,"BM_ALLOC_ROOT")

	alls=leaf+sap_ab+sap_be+fruit+root

	biomass=ncvar_get(new,"TOTAL_M")
	st=ncvar_get(new,"ST_MORTALITY")
	dd=ncvar_get(new,"TOTAL_BM_LITTER")

	leaf_turn=ncvar_get(new,"LEAF_TURN")
	sap_turn=ncvar_get(new,"SAP_AB_TURN")
	root_turn=ncvar_get(new,"ROOT_TURN")
	fruit_turn=ncvar_get(new,"FRUIT_TURN")

	all_turn=leaf_turn+sap_turn+root_turn+fruit_turn

	st=ncvar_get(new,"ST_MORTALITY")
	
	leaf=ncvar_get(new,'LEAF_BM_LITTER')
	sap=ncvar_get(new,'SAP_AB_BM_LITTER')
	heart=ncvar_get(new,'HEART_AB_BM_LITTER')
	fruit=ncvar_get(new,'FRUIT_BM_LITTER')
	sap.be=ncvar_get(new,'SAP_BE_BM_LITTER')
	heart.be=ncvar_get(new,'HEART_BE_BM_LITTER')
	carbre=ncvar_get(new,'RESERVE_BM_LITTER')
	root=ncvar_get(new,'ROOT_BM_LITTER')
	dd=leaf+sap+heart+fruit+sap.be+heart.be+root+carbre
	
	#dd=ncvar_get(new,"TOTAL_BM_LITTER")

	biomass=ncvar_get(new,"TOTAL_M")
	

	real_mortality=rep(0,2920)
	for (i in 1:2920) {
		if (dd[i]==0) {
			if (st[i]!=0) {
				real_mortality[i]=st[i]
			} else {
				real_mortality[i]=0
			}
		}   else {
			if (dd[i]>=st[i]) {
				real_mortality[i]=dd[i]
			} else if (dd[i]<st[i]) {
				real_mortality[i]=st[i]
			}
		}
	}
	expe_data[,ex,1]=alls
	expe_data[,ex,2]=all_turn
	expe_data[,ex,3]=st
	expe_data[,ex,4]=dd
	expe_data[,ex,5]=real_mortality
	
	

	env_mortality=real_mortality-st
	expe_data[,ex,6]=env_mortality
	expe_data[,ex,7]=biomass
	
	for (year in 1:8) {
		expe_sum[year,ex,1]=sum(alls[((year-1)*365+1):(year*365)])
		expe_sum[year,ex,2]=sum(all_turn[((year-1)*365+1):(year*365)])
		expe_sum[year,ex,3]=sum(st[((year-1)*365+1):(year*365)])
		expe_sum[year,ex,4]=sum(dd[((year-1)*365+1):(year*365)])
		expe_sum[year,ex,5]=sum(real_mortality[((year-1)*365+1):(year*365)])
		expe_sum[year,ex,6]=sum(env_mortality[((year-1)*365+1):(year*365)])
		expe_sum[year,ex,7]=biomass[(year*365)]-biomass[((year-1)*365+1)]
	}
}


expe_sum=expe_sum/100


load('biomass_mortality.Rdata')

tiff("stacked delta biomass 11June.tiff",width=800,height=800)
par(mfrow=c(2,1))
par(mar=c(5,5,3,3))
ex=1
plot(0,0,col="white",xlim=c(0.5,8.5),xlab="",xaxt="n",ylab=expression('△biomass'~(MgC~ha^-1~yr^-1)),cex.lab=2,cex.axis=2,main="CTL",cex.main=3,
ylim=c(-(max(expe_sum[,,5]+expe_sum[,,2])),max(expe_sum[,,1])))

for (i in 1:8) {
	rect(i-0.3,0,i+0.3,expe_sum[i,ex,1],col="#7fbf7b",border=NA)
}
for (i in 1:8) {
	rect(i-0.3,-expe_sum[i,ex,2],i+0.3,0,col="#ef8a62",border=NA)
}
for (i in 1:8) {
	rect(i-0.3,-(expe_sum[i,ex,2]+expe_sum[i,ex,3]),i+0.3,-expe_sum[i,ex,2],col="#af8dc3",border=NA)
}
for (i in 1:8) {
	rect(i-0.3,-(expe_sum[i,ex,2]+expe_sum[i,ex,3]+expe_sum[i,ex,6]),i+0.3,-(expe_sum[i,ex,2]+expe_sum[i,ex,3]),
	col="#b2182b",border=NA)
}
abline(h=0,lty=2)
lines(expe_sum[,ex,1]-expe_sum[,ex,2]-expe_sum[,ex,5],lwd=4)
axis(1, at=seq(-1,by=1,length.out=12), 
						labels=seq(1999,by=1,length.out=12),cex.axis=2,lwd=2)
axis(2,at=seq(-30,20,by=10),labels=FALSE,lwd=2)
axis(3, at=seq(-1,by=1,length.out=12), labels=FALSE,lwd=2,tck=0)
axis(4,at=seq(-30,20,by=10),labels=FALSE,lwd=2,tck=0)

points(-biomass_mortality[1:8,1],pch=19,lwd=2,cex=2)
for (yy in 1:8) {
	top=biomass_mortality[16+yy,1]
	med=biomass_mortality[yy,1]
	down=2*med-top
	lines(c(yy,yy),c(-top,-down),lwd=2)
}

ex=2
par(mar=c(5,5,3,3))
plot(0,0,col="white",xlim=c(0.5,8.5),xlab="",xaxt="n",ylab=expression('△biomass'~(MgC~ha^-1~yr^-1)),cex.lab=2,cex.axis=2,,main="TFE",cex.main=3,
ylim=c(-(max(expe_sum[,,5]+expe_sum[,,2])),max(expe_sum[,,1])))
for (i in 1:8) {
	rect(i-0.3,0,i+0.3,expe_sum[i,ex,1],col="#7fbf7b",border=NA)
}
for (i in 1:8) {
	rect(i-0.3,-expe_sum[i,ex,2],i+0.3,0,col="#ef8a62",border=NA)
}
for (i in 1:8) {
	rect(i-0.3,-(expe_sum[i,ex,2]+expe_sum[i,ex,3]),i+0.3,-expe_sum[i,ex,2],col="#af8dc3",border=NA)
}
for (i in 1:8) {
	rect(i-0.3,-(expe_sum[i,ex,2]+expe_sum[i,ex,3]+expe_sum[i,ex,6]),i+0.3,-(expe_sum[i,ex,2]+expe_sum[i,ex,3]),
	col="#b2182b",border=NA)
}
abline(h=0,lty=2)
lines(expe_sum[,ex,1]-expe_sum[,ex,2]-expe_sum[,ex,5],lwd=4)
axis(1, at=seq(-1,by=1,length.out=12), 
						labels=seq(1999,by=1,length.out=12),cex.axis=2,lwd=2)
axis(2,at=seq(-30,20,by=10),labels=FALSE,lwd=2)
axis(3, at=seq(-1,by=1,length.out=12), labels=FALSE,lwd=2,tck=0)
axis(4,at=seq(-30,20,by=10),labels=FALSE,lwd=2,tck=0)

points(-biomass_mortality[9:16,1],pch=19,lwd=2,cex=2)
for (yy in 1:8) {
	top=biomass_mortality[24+yy,1]
	med=biomass_mortality[yy+8,1]
	down=2*med-top
	lines(c(yy,yy),c(-top,-down),lwd=2)
}

legend("bottomright",legend=c("growth","turnover","self-thinning","climate","△biomass"),
text.col=c("#7fbf7b","#ef8a62","#af8dc3","#b2182b","black"),bty="n",cex=2)
legend("bottomleft",expression("biomass mortality \nda Costa et al 2010"),pch=19,bty="n",cex=2)


dev.off()
