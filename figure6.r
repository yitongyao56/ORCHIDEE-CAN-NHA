source('default.R')
modelname=c("11NOV1129usdav130223kill17gmax1600p50leaf",
"11NOV1129usdav130223kill17gmax1600p50leaf",
"11NOV1129usdav130223kill17gmax1600",
"11NOV0606gpsi1kf20copygmax1400day15",
"emilieTEST0427MR001minus17keepDD0430",
"11NOVdrain04v41129","11NOV1129")
all_soil=c("moistc_2","phi_soil_2")
soildata=rep(NA,96*5*2*2*12)
dim(soildata)=c(96,5,2,2,12)
for (num in 4:4) {
	model=modelname[num]
	for (ex in 1:2)  {
		if (num==2) {
			new=nc_open(paste0("CAXclose",expe[ex],"newKsap5startDD",model,"sechiba.nc"))
		} else if (num==1) {
			new=nc_open(paste0("CAXclose2",expe[ex],"newKsap5startDD",model,"sechiba.nc"))
		} else {
			new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"sechiba.nc"))
		}	
		for (vars in 1:2) {
			bles=ncvar_get(new,all_soil[vars])
			for (layer in 1:12) {
				mo_var=myfunction(bles[layer,])
				soildata[,num,ex,vars,layer]=mo_var
			}
		}
	}
}

# par(mfrow=c(3,2))
library(RColorBrewer)

library(shape) # for Arrows
solay=ncvar_get(new,"solay")


tiff(paste0("smc heatmap six plots nopoints 11June.tiff"),height=500,width=1200)
par(fig=c(0.1,0.9,0.05,0.15))
#plot(1:10,1:10,col="white",axes=FALSE)
par(mar=c(2,2,2,3))
colos=brewer.pal(10, "Spectral")
image.plot(legend.only=TRUE,zlim=c(0.05,0.3),col=colos,horizontal=TRUE,
breaks=c(0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3),legend.mar=1,
legend.width=8,legend.shrink=1,axis.args=list(cex.axis=2)) #legend.args=list(cex=2,text=expression(SMC (m^"3"/m^"3"))))
mtext(expression(SMC~(m^"3"~m^"-3")),cex=2,side=1,adj=0.5,padj=-1)
depth=c(0,0.5,1,2,4)
tolayer=12-c(0.5,9,10,11)
for (model in 4:4) {
	for (ex in 1:2) {
		par(fig=c(0.04+(ex-1)*0.43,0.04+ex*0.43,0.2,0.8),new=TRUE)
		par(mar=c(2,5,2,2))
		df=soildata[,model,ex,1,]
		colos=brewer.pal(10, "Spectral")
		image(df[,12:1],las=1,zlim=c(0.05,0.3),
		breaks=c(0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3),
		col=colos,xaxt="n",ylab="",cex.axis=2,
		cex.main=2,yaxt="n",legend.cex=2,cex.lab=2,legend.width=2,
		main=expe[ex],legend.shrink=1,ylim=c(0,1),
		axis.args=list(cex.axis=2),yaxs="i")
		axis(1,at=seq(-1/8,by=1/8,length.out=10),labels=c(" "," ",seq(from=2002,to=2009,by=1)),cex.axis=2,lwd=2)
		axis(2,at=seq(-1/12,by=1/12,length.out=13),
		labels=c(' ',round(solay,3)[12:1]),cex.axis=2,lwd=2,las=1)
		Arrows(x0 = 1/8, y0 = 1.5, x1 = 1/8, y1 = 1.03,xpd=TRUE,lwd=2,arr.width=0.15)	
		axis(3,at=seq(-1/8,by=1/8,length.out=10),,cex.axis=2,lwd=2,tck=0,labels=FALSE)
		axis(4,at=seq(-1/12,by=1/12,length.out=13),tck=0,lwd=2,labels=FALSE)
		if (ex==1) {
		mtext('Depth (m)',side=2,cex=2,adj=0.5,line=6)
		}
		if (ex==1) {	
			for (aa in 1:4) {
				ev_text=paste0("obs_smc_ctl=soilwater$mean_ctl",depth[aa],"m")
				eval(parse(text=ev_text))	
				obs_smc_ctl=obs_smc_ctl/100
				
				tt=time_f(soilwater$X)/96
				obs_cut=colos[as.numeric(cut(obs_smc_ctl,
				breaks=c(0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3)))]
				#points(tt,rep(tolayer[aa]/12,length(tt)),pch=21,bg=obs_cut,col='gray',cex=1.5)
			
			}
		} else {
			for (aa in 1:4) {
				ev_text=paste0("obs_smc_tfe=soilwater$mean_tfe",depth[aa],"m")
				eval(parse(text=ev_text))	
				obs_smc_tfe=obs_smc_tfe/100
				tt=time_f(soilwater$X)/96
				obs_cut=colos[as.numeric(cut(obs_smc_tfe,
				breaks=c(0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3)))]
				#points(tt,rep(tolayer[aa]/12,length(tt)),pch=21,bg=obs_cut,col='gray',cex=1.5)
			}
		}
		
		# depth 
	}
}
dev.off()

