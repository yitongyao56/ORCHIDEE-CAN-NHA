source('default.R')

vl_name=c("PSLeaf","PSStem","PSRoot","stemPLC")

modelname="11NOV0606gpsi1kf20copygmax1400day15"
vl_month_data=rep(NA,96*2*length(vl_name)*20)
dim(vl_month_data)=c(96,2,length(vl_name),20)
vl_daily_data=rep(NA,2920*2*length(vl_name)*20)
dim(vl_daily_data)=c(2920,2,length(vl_name),20)
for (ex in 1:2) {

	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname,"sechiba.nc"))   
	for (num in 1:length(vl_name)) {
		for (co in 1:20) {
			vardata=ncvar_get(new,paste0(vl_name[num],"_",sprintf("%03d",co)))
			mon=myfunction(vardata)
			vl_month_data[,ex,num,co]=mon
			vl_daily_data[,ex,num,co]=vardata
		}
		
	}
}

heightdata=rep(NA,96*2*20)
dim(heightdata)=c(96,2,20)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname,"stomate.nc"))
	for (co in 1:20) {
		vardata=ncvar_get(new,paste0("CCH_",sprintf("%03d",co)))
		mon=myfunction(vardata)
		heightdata[,ex,co]=mon
	}
}

# co 5 10 15 
ileaf=1
istem=2
iroot=3
co=5
wet=53  # 53
dry=59# 59

library(RColorBrewer)
colos=brewer.pal(9,"Spectral")
#var_label=c("psi_leaf","psi_stem","psi_root")
var_label=c("Ψleaf","Ψstem","Ψroot","Ψ gradient")
ylabs=c('','Cohorts')
tiff("figure 8 phi organs ncirc 9June.tiff",width=1300,height=600)
par(mfrow=c(2,4))
for (ex in 1:2) {
	for (vars in 1:3) {
		#par(fig=c(0+(vars-1)*0.33,vars*0.33,0.99-(0.49*ex),0.99-0.49*(ex-1)),new=TRUE)
		par(mar=c(5,5,3,6))
		df=vl_month_data[,ex,vars,]
		df[df<(-2.5)]=-2.5
		if (vars==1) {
			ylabs='Cohorts'
		} else {
			ylabs=''
		}
		image.plot(df[,20:1],las=1,zlim=c(-2.5,0),col=colos,xaxt="n",ylab=ylabs,cex.axis=2,cex.lab=2,
		cex.main=3,yaxt="n",legend.cex=2,legend.mar=11,cex.lab=2,legend.width=2,legend.shrink=1,
		axis.args=list(cex.axis=2),main=var_label[vars])
		#par(fig=c(0+(vars-1)*0.33,vars*0.33,0.99-(0.49*ex),0.99-0.49*(ex-1)),new=TRUE)
		
		Arrows(x0 = 1/8, y0 = 1.5, x1 = 1/8, y1 = 1.03,xpd=TRUE,lwd=2,arr.width=0.15)
		axis(1,at=seq(-1/8,by=1/8,length.out=11),labels=seq(2000,by=1,length.out=11),cex.axis=2,lwd=2)
		axis(2,at=seq(-1/10,by=1/10,length.out=11),labels=c(' ',seq(20,by=-2,length.out=10)),cex.axis=2,lwd=2)	
		axis(3,at=seq(-1/8,by=1/8,length.out=11),labels=FALSE,cex.axis=2,lwd=2,tck=0)
		axis(4,at=seq(-1/10,by=1/10,length.out=11),labels=FALSE,cex.axis=2,lwd=2,tck=0)	
		if (ex==1 & vars==1) {
			mtext("CTL",cex=2,side=3,line=0.8,at=-0.1)
		}
		if (ex==2 & vars==1) {
			mtext("TFE",cex=2,side=3,line=0.8,at=-0.1)
		}
		mtext("MPa",cex=2,side=3,line=1,at=1.1)
	}
	par(mar=c(5,5,3,6))
	pchs=c(5,19)
	plot(1:10,1:10,pch=19,col="white",xlab="Ψ (MPa)",ylab="",yaxt='n',
	cex.lab=2,cex.axis=2,cex=2,xlim=c(-2.5,0),ylim=c(0,40),main=var_label[4],cex.main=2)
	axis(4,at=seq(-10,50,by=10),cex.axis=2,lwd=2)
	axis(2,at=seq(-10,50,by=10),labels=FALSE,cex.axis=2,lwd=2,tck=0)
	axis(1,seq(-3,0.5,by=0.5),labels=FALSE,lwd=2)
	axis(3,seq(-3,0.5,by=0.5),labels=FALSE,lwd=2,tck=0)
	
	mtext("Height (m)", side=4, line=3.5,cex=1.5)

	aa=0
	for (co in c(5,10)) {
		wet_psi_pro=c(vl_month_data[wet,ex,ileaf,co],vl_month_data[wet,ex,istem,co],vl_month_data[wet,ex,iroot,co])
		wet_height_pro=c(heightdata[wet,ex,co],heightdata[wet,ex,co]/2,0)

		dry_psi_pro=c(vl_month_data[dry,ex,ileaf,co],vl_month_data[dry,ex,istem,co],vl_month_data[dry,ex,iroot,co])
		dry_height_pro=c(heightdata[dry,ex,co],heightdata[dry,ex,co]/2,0)
		
		aa=aa+1
		points(wet_psi_pro,wet_height_pro,pch=pchs[aa],cex=2)	
		lines(wet_psi_pro,wet_height_pro,lty=1,lwd=2)
		points(dry_psi_pro,dry_height_pro,col="red",pch=pchs[aa],cex=2)
		lines(dry_psi_pro,dry_height_pro,lty=1,col="red",lwd=2)
		if (ex==2) {
			labs=c('L','S','R')
			if (co==5) {			
				for (jj in 1:3) {
					text(dry_psi_pro[jj]-0.2,dry_height_pro[jj],labs[jj],cex=2)
				}
			} else {
				labs=c('L','S','R')
				for (jj in 1:3) {
					text(dry_psi_pro[jj]+0.2,dry_height_pro[jj],labs[jj],cex=2)
				}
			}
		}
	}
	if (ex==2) {
		legend("topright",c("wet season","dry season"),text.col=c("black","red"),cex=2,bty="n")
		legend("right",c("#5","#10"),bty="n",pch=c(5,19),cex=2)
	}
	
}

dev.off()