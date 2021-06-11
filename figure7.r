
source('default.R')
modelname=c("falseROOTfalseHYDRO","trueROOTfalseHYDRO",
"11NOV0606gpsi1kf20copygmax1400day15") 

carbon_vars=c("GPP","LAI","IND","TOTAL_M")
carbon_data=rep(NA,2920*3*2*4)
dim(carbon_data)=c(2920,3,2,4)

carbon_daily=rep(NA,96*3*2*4)
dim(carbon_daily)=c(96,3,2,4)
for (num in 3:3) {
	model=modelname[num]
	for (ex in 1:2) {
		new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",model,"stomate.nc"))
		for (vars in 1:4) {     # 19maint
			vard=ncvar_get(new,carbon_vars[vars])
			carbon_data[,num,ex,vars]=vard
			mon=myfunction(vard)
			carbon_daily[,num,ex,vars]=mon
		}
	}
}


library(R.matlab)
cax=readMat("caxGPPdata.mat")
caxGPP=cax$caxGPP

CAXgppfunction=function(stom_var) {
   month_var=c()
   for (year in 1:3) {
        year_series=stom_var[((year-1)*365+1):(year*365)]
        for (month in 1:12) {
           month_series= year_series[(cum_monthdays[month]-monthdays[month]+1):cum_monthdays[month]]
           month_var=c(month_var,mean(month_series))
        }
   }
   return(month_var)
}
ctlGPP=read.csv("CTL_GPP.csv",header=F)
tfeGPP=read.csv("TFE_GPP.csv",header=F)

fluxGPP=CAXgppfunction(caxGPP)


monthdays=c(31,28,31,30,31,30,31,31,30,31,30,31)
cum_monthdays=cumsum(monthdays)
# lai gpp biomass Re
ctlgppfunction=function(stom_var) {
   month_var_ctl=c()
   for (year in 1:2) {
		find_row=(stom_var[,1]>(year*365)) & (stom_var[,1]<=((year*365)+365))
		varvar=stom_var[find_row,]
        year_series=varvar[,1]-(year*365)
        for (month in 1:12) {
           month_series= (year_series>=(cum_monthdays[month]-monthdays[month]+1)) & (year_series<=cum_monthdays[month])
           month_var_ctl=c(month_var_ctl,mean(varvar[month_series,2],na.rm=T))
        }
   }
   return(month_var_ctl)
}
monthdays=c(31,28,31,30,31,30,31,31,30,31,30,31)
cum_monthdays=cumsum(monthdays)
# lai gpp biomass Re
tfegppfunction=function(stom_var) {
   month_var_ctl=c()
   for (year in 0:2) {
		find_row=(stom_var[,1]>(year*365)) & (stom_var[,1]<=((year*365)+365))
		varvar=stom_var[find_row,]
        year_series=varvar[,1]-(year*365)
        for (month in 1:12) {
           month_series= (year_series>=(cum_monthdays[month]-monthdays[month]+1)) & (year_series<=cum_monthdays[month])
           month_var_ctl=c(month_var_ctl,mean(varvar[month_series,2],na.rm=T))
        }
   }
   return(month_var_ctl)
}

obs_ctlGPP=ctlgppfunction(ctlGPP)
obs_tfeGPP=tfegppfunction(tfeGPP)


tiff("figure 7 GPP 11June.tiff",width=1000,height=400)
par(mfrow=c(1,2))
par(mar=c(5,6,3,3))
yyy=c(c(obs_tfeGPP[1:12],obs_ctlGPP),carbon_daily[1:36,3,1,1],fluxGPP)
plot(carbon_daily[1:36,3,1,1],type="l",ylim=c(1,13),lwd=2,
ylab=expression(GPP~(gC~m^-2~d^-1)),xlab="",xaxt="n",cex.axis=2,cex.lab=2,main="CTL",cex.main=2)
rect(7,0,11,15,col="#d9d9d9",border=NA)
rect(19,0,23,15,col="#d9d9d9",border=NA)
rect(31,0,35,15,col="#d9d9d9",border=NA)
lines(carbon_daily[1:36,3,1,1],lwd=2)
lines(c(obs_tfeGPP[1:12],obs_ctlGPP),lwd=2,lty=2)
lines(7:36,fluxGPP[7:36],lwd=2,lty=3)
Arrows(x0 = 13, y0 = -3, x1 = 13, y1 = -1,xpd=TRUE,lwd=2,arr.width=0.15)
axis(1, at=seq(-11,by=12,length.out=5), 
					labels=seq(2000,by=1,length.out=5),cex.axis=2,lwd=2)
legend("bottomleft",legend=c("ORCHIDEE-CAN-NHA","SPA model","flux observation (only CTL)"),
lty=c(1,2,3),cex=2,bty="n")	
axis(2,at=seq(0,14,by=2),lwd=2,labels=FALSE)
axis(3,at=seq(1,by=12,length.out=4),cex.axis=2,lwd=2,tck=0,labels=FALSE)
axis(4,at=seq(0,14,by=2),tck=0,lwd=2,labels=FALSE)

par(mar=c(5,6,3,3))
yyy=c(obs_tfeGPP,carbon_daily[1:36,3,2,1])
plot(carbon_daily[1:36,3,2,1],type="l",ylim=c(1,13),lwd=2,
ylab=expression(GPP~(gC~m^-2~d^-1)),xlab="",xaxt="n",cex.axis=2,cex.lab=2,main="TFE",cex.main=2)
rect(7,0,11,15,col="#d9d9d9",border=NA)
rect(19,0,23,15,col="#d9d9d9",border=NA)
rect(31,0,35,15,col="#d9d9d9",border=NA)
lines(carbon_daily[1:36,3,2,1],lwd=2)
lines(obs_tfeGPP,lwd=2,lty=2)
#legend("bottomleft",legend=c("ORCHIDEE-CAN-NHA","SPA model"),
#lty=c(1,2),cex=2,bty="n")
Arrows(x0 = 13, y0 = -3, x1 = 13, y1 = -1,xpd=TRUE,lwd=2,arr.width=0.15)
axis(1, at=seq(-11,by=12,length.out=5), 
						labels=seq(2000,by=1,length.out=5),cex.axis=2)
axis(2,at=seq(0,14,by=2),lwd=2,labels=FALSE)
axis(3,at=seq(1,by=12,length.out=4),cex.axis=2,lwd=2,tck=0,labels=FALSE)
axis(4,at=seq(0,14,by=2),tck=0,lwd=2,labels=FALSE)

dev.off()