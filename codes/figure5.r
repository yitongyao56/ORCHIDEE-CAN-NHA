source('default.R')
time_f=function(tt) {
	time_var=c()
	for (year in 2001:2004) {
		aa=tt[tt>=((year-2001)*365+1) & tt<=((year-2001)*365+365)]
		aa=aa-((year-2001)*365)
		time_var=c(time_var,((aa/365)*12)+((year-2001)*12))
	}	
	return(time_var)
}
smc_name=c("0m","0.5m","1m","2m","4m")
soilwater=read.csv("ESECAFLOR_soil_water_5m_summary_2000_2004.csv",header=T)
soilwater=soilwater[1:65,]
model_depth=c("005","05","1","2","4")
depth=c(0,0.5,1,2,4)
soilwater=soilwater[(soilwater$X)>=0 & !is.na(soilwater$X),]

model="Ksap5"
all_soil=c("moistc_2")
soildata=rep(NA,96*2*2*12)
dim(soildata)=c(96,2,2,12)

for (ex in 1:2)  {
	new=nc_open(paste0("CAX",expe[ex],
	"newKsap5startDD11NOV0606gpsi1kf20copygmax1400day15sechiba.nc"))
	vars=1
	bles=ncvar_get(new,all_soil[vars])
	for (layer in 1:12) {
		mo_var=myfunction(bles[layer,])
		soildata[,ex,1,layer]=mo_var
	}
}



tiff("figure 5 smc 11June.tiff",width=1200,height=500)
par(mfrow=c(2,4))
colos=c("black")
aa=0
ex=1

for (layer in c(1,9,10,11)) {
	aa=aa+1
	vars=2
	par(mar=c(5,7,4,2))
	if (vars==2 & layer==1) {
		plot(soildata[1:48,ex,1,layer],type="l",col="white",xaxt="n",xlab="",ylab=expression(SMC~(m^3~m^-3)),
		cex.lab=2,cex.axis=2,main=smc_name[aa],cex.main=2,ylim=c(0.00,0.40),yaxs="i")
	}
	if (vars==2 & layer!=1) {
		plot(soildata[1:48,ex,1,layer],type="l",col="white",xaxt="n",xlab="",ylab="",
		cex.lab=2,cex.axis=2,main=smc_name[aa],cex.main=2,ylim=c(0.00,0.40),yaxs="i")
	}
	axis(2, at=seq(-0.1,0.5,by=0.1), labels = FALSE,cex.axis=2, lwd = 2)
	#axis(1,at=seq(0,6,by=1), labels = FALSE,lwd = 2)
	axis(1, at=seq(-11,by=12,length.out=9), 
			labels=seq(2000,by=1,length.out=9),cex.axis=2,lwd=2)	
	axis(3,at=seq(-11,by=12,length.out=9),cex.axis=2,lwd=2,tck=0,labels=FALSE)
	axis(4,at=seq(-0.1,0.5,by=0.1),tck=0,lwd=2,labels=FALSE)	
	
	
	lines(soildata[1:48,ex,1,layer],lwd=2)
	rect(7,-1,11,15,col="#d9d9d9",border=NA)
	rect(19,-1,23,15,col="#d9d9d9",border=NA)
	rect(31,-1,35,15,col="#d9d9d9",border=NA)
	rect(43,-1,47,15,col="#d9d9d9",border=NA)
	#for (model in c(1,2,3,4)) {
		lines(soildata[1:48,ex,1,layer],lwd=2,col=colos[1])
		#lines(soildata[1:48,ex,2,layer],lwd=2,col="blue")
	#}
	# add observation 
	ev_text=paste0("obs_smc_ctl=soilwater$mean_ctl",depth[aa],"m")
	eval(parse(text=ev_text))	
	ev_text=paste0("std_smc_ctl=soilwater$std_ctl",depth[aa],"m")
	eval(parse(text=ev_text))
	obs_smc_ctl=obs_smc_ctl/100
	std_smc_ctl=std_smc_ctl/100
	if (vars==2) {
		tt=time_f(soilwater$X)
		points(tt,obs_smc_ctl,pch=1,col="black")
		errorbar(tt,obs_smc_ctl,0,std_smc_ctl,bar.col="black",add=TRUE)
	}
	Arrows(x0 = 13, y0 = 0.06, x1 = 13, y1 = 0.02,xpd=TRUE,lwd=2,arr.width=0.15)
	
	axis(1, at=seq(1,by=12,length.out=5), 
			labels=seq(2001,by=1,length.out=5),cex.axis=2)	
	if (aa==1) {
		mtext("CTL",cex=2,side=3,line=1.5,at=-8)
	}
}

aa=0
ex=2
for (layer in c(1,9,10,11)) {
	aa=aa+1
	vars=2
	par(mar=c(5,7,4,2))
	# expression(paste0("SMC (m"^"3""/m"^"3"")"))  "SMC (m"^"3"*"/m"^"3"*")"
	if (vars==2 & layer==1) {
		plot(soildata[1:48,ex,1,layer],type="l",col="white",xaxt="n",xlab="",ylab=expression(SMC~(m^3~m^-3)),
		cex.lab=2,cex.axis=2,main=smc_name[aa],cex.main=2,ylim=c(0.00,0.40),yaxs="i")
	}
	if (vars==2 & layer!=1) {
		plot(soildata[1:48,ex,1,layer],type="l",col="white",xaxt="n",xlab="",ylab="",
		cex.lab=2,cex.axis=2,main=smc_name[aa],cex.main=2,ylim=c(0.00,0.40),yaxs="i")
	}
	axis(2, at=seq(-0.1,0.5,by=0.1), labels = FALSE,cex.axis=2, lwd = 2)
	#axis(1,at=seq(0,6,by=1), labels = FALSE,lwd = 2)
	axis(1, at=seq(-11,by=12,length.out=9), 
			labels=seq(2000,by=1,length.out=9),cex.axis=2,lwd=2)	
	axis(3,at=seq(-11,by=12,length.out=9),cex.axis=2,lwd=2,tck=0,labels=FALSE)
	axis(4,at=seq(-0.1,0.5,by=0.1),tck=0,lwd=2,labels=FALSE)	
	
	
	lines(soildata[1:48,ex,1,layer],lwd=2)
	rect(7,-1,11,15,col="#d9d9d9",border=NA)
	rect(19,-1,23,15,col="#d9d9d9",border=NA)
	rect(31,-1,35,15,col="#d9d9d9",border=NA)
	rect(43,-1,47,15,col="#d9d9d9",border=NA)
	#for (model in c(1,2,3,4)) {
		lines(soildata[1:48,ex,1,layer],lwd=2,col=colos[1])
		#(soildata[1:48,ex,2,layer],lwd=2,col="blue")
	#}
	# add observation 
	ev_text=paste0("obs_smc_tfe=soilwater$mean_tfe",depth[aa],"m")
	eval(parse(text=ev_text))	
	ev_text=paste0("std_smc_tfe=soilwater$std_tfe",depth[aa],"m")
	eval(parse(text=ev_text))
	obs_smc_tfe=obs_smc_tfe/100
	std_smc_tfe=std_smc_tfe/100
	if (vars==2) {
		tt=time_f(soilwater$X)
		points(tt,obs_smc_tfe,pch=1,col="black")
		errorbar(tt,obs_smc_tfe,0,std_smc_tfe,bar.col="black",add=TRUE)
	}
	Arrows(x0 = 13, y0 = 0.06, x1 = 13, y1 = 0.02,xpd=TRUE,lwd=2,arr.width=0.15)
	axis(1, at=seq(1,by=12,length.out=5), 
			labels=seq(2001,by=1,length.out=5),cex.axis=2)	
	if (aa==1) {
		mtext("TFE",cex=2,side=3,line=1.5,at=-8)
	}
}

dev.off()
