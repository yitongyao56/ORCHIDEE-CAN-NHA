source('default.R')
sap=read.csv("sap_flow_summary_2002_2003.csv",header=T)
sap=sap[,1:4]
names(sap)=c("Date","Period","CTL","TFE")

monthdays=c(31,28,31,30,31,30,31,31,30,31,30,31)
cum_monthdays=cumsum(monthdays)
# lai gpp biomass Re
sapfunction=function(stom_var) {
   month_var_ctl=data.frame()
   month_var_tfe=data.frame()
   #names(month_var_ctl)="ctl"
   #names(month_var_tfe)="tfe"
   for (year in 1:2) {
		find_row=(stom_var[,2]>(year*365)) & (stom_var[,2]<=((year*365)+365))
		varvar=stom_var[find_row,]
        year_series=varvar[,2]-(year*365)
        for (month in 1:12) {
           month_series= (year_series>=(cum_monthdays[month]-monthdays[month]+1)) & (year_series<=cum_monthdays[month])
           month_var_ctl=rbind(month_var_ctl,mean(varvar[month_series,3],na.rm=T))
		   month_var_tfe=rbind(month_var_tfe,mean(varvar[month_series,4],na.rm=T))
        }
   }
   month_var=cbind(month_var_ctl,month_var_tfe)
   return(month_var)
}
month_sap=sapfunction(sap)

expe=c("CTL","TFE") # 11NOV1129usdav13
modelname=c("","11NOV0606gpsi1kf20copygmax1400day15") 
supplyname=c("transpir_supply","transpir_tot_tzj")

Edata=rep(NA,2*2*2*96)
dim(Edata)=c(2,2,2,96)
for (model in 1:2) {
	for (ex in 1:2) {
		new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname[model],"sechiba.nc"))
		transpir=ncvar_get(new,"transpir")
		supply=ncvar_get(new,"transpir_supply")
		mon_trans=myfunction(transpir)
		mon_supply=myfunction(supply)
		Edata[model,ex,1,]=mon_trans
		Edata[model,ex,2,]=mon_supply
		#nc_close(new)
	}
}


ann_tmp=rep(NA,365*8)
dim(ann_tmp)=c(365*8,1)
max_tmp=rep(NA,365*8)
dim(max_tmp)=c(365*8,1)
ann_prec=rep(NA,365*8)
dim(ann_prec)=c(365*8,1)
ann_swdown=rep(NA,365*8)
dim(ann_swdown)=c(365*8,1)

for (year in 2001:2008) {
	force=nc_open(paste0("./caxiuana/caxiuana_",year,".nc"))
	tair=ncvar_get(force,"Tair")
	rain=ncvar_get(force,"Rainf")
	swdown=ncvar_get(force,"SWdown")
	rain=rain*1800 # ->0.5h
	for (day in 1:365) {
		vpd=mean(tair[((day-1)*48+1):(day*48)],na.rm=T)
		ann_tmp[(year-2001)*365+day,1]=vpd-273.15
		vpd=max(tair[((day-1)*48+1):(day*48)],na.rm=T)
		max_tmp[(year-2001)*365+day,1]=vpd-273.15
		vpd=sum(rain[((day-1)*48+1):(day*48)],na.rm=T)
		ann_prec[(year-2001)*365+day,1]=vpd
		vpd=mean(swdown[((day-1)*48+1):(day*48)],na.rm=T)
		ann_swdown[(year-2001)*365+day,1]=vpd
	}	
}


modelname=c("","11NOV0606gpsi1kf20copygmax1400day15")  
deficit_model=rep(NA,96*2*2)
dim(deficit_model)=c(96,2,2)
for (model in 1:2) {
	for (ex in 1:2) {
		annual_prec=ann_prec
		month_prec=myfunction(annual_prec)
		new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD",modelname[model],"sechiba.nc"))
		evap=ncvar_get(new,"evap")
		month_evap=myfunction(evap)
		if (ex==2) {
			month_prec[13:96]=month_prec[13:96]/2
		}
		deficit=month_prec-month_evap
		deficit_model[,model,ex]=deficit
	}
}


for (model in 1:2) {
	for (ex in 1:2) {
		obs=month_sap[,ex]
		preds=Edata[model,ex,1,13:36]
		rows=which(!is.na(obs))
		print(cor(obs[rows],preds[rows]))
	}
}

library(Metrics)
for (model in 1:2) {
	for (ex in 1:2) {
		obs=month_sap[,ex]
		preds=Edata[model,ex,1,13:36]
		rows=which(!is.na(obs))
		print(rmse(obs[rows],preds[rows]))
	}
}


Rcor=c(0.77,0.51,0.80,0.51)
aa=0
tiff("figure 3 11June.tiff",width=700,height=700)
model_label=c("ORCHIDEE-CAN-RS","ORCHIDEE-CAN-NHA")

deficit_model[deficit_model>0]=0
colos=brewer.pal(9,"YlOrRd")
colos=colos[9:1]
point.cols=colos[as.numeric(cut(deficit_model,breaks=c(-3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0)))]
breakss=c(-3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0)
# model 1 0.05-0.5 model 2 0.5 0.95 
# model 1 0.5 0.95 model 2 0.05-0.5

par(mfrow=c(2,2))
for (model in 1:2) {
	if (model==2) {
		xlabs=expression(Observed~sap~flow~(mm~d^-1))
	} else {
		xlabs=""
	}
	for (ex in 1:2) {
		aa=aa+1
		par(fig=c((ex-1)*0.45+0.02,ex*0.45+0.02,1-(model*0.45+0.05),1-((model-1)*0.45+0.06)),new=TRUE)
		par(mar=c(4,5,3,2))
		yy=c(month_sap[,ex],Edata[model,ex,1,13:36])
		defo=deficit_model[13:36,model,ex]
		defo[defo<(-3)]=-2.9
		point.cols=colos[as.numeric(cut(defo,breaks=c(-3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0)))]
		plot(month_sap[,ex],Edata[model,ex,1,13:36],
		bg=point.cols,cex=2,pch=21,ylim=c(1,4),
		xlim=c(1,4),main=paste0(expe[ex]),cex.main=2,ylab=paste0(model_label[model]," (mm/d)"),
		xlab=xlabs,cex.lab=2,cex.axis=2)	
		##
		axis(2, at=seq(-1,6,by=1), labels = FALSE,cex.axis=1.25, lwd = 2)
		axis(1,at=seq(-1,6,by=1), labels = FALSE,lwd = 2)
		axis(3,at=seq(-1,6,by=1),tck=0,lwd=2,labels=FALSE)
		axis(4,at=seq(-1,6,by=1),tck=0,lwd=2,labels=FALSE)	
		##
		lmm=lm(Edata[model,ex,1,13:36] ~ month_sap[,ex])
		xx=seq(from=min(yy,na.rm=T),to=max(yy,na.rm=T),length=20)
		newy=xx*lmm$coefficients[2]+lmm$coefficients[1]
		lines(xx,newy,col="red",lty=2,lwd=2)
		abline(a=0,b=1,lty=2,col="black",lwd=2)
		legend("bottomright",legend=paste0("R=",Rcor[aa]),cex=2.5,bty="n")
	}
	if (ex==2) {
		# image.plot(legend.only=TRUE,zlim=c(-4.5,0),col=colos,breaks=c(-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0),
		# lab.breaks=c(-3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0),legend.shrink = 1, legend.width = 1.2)
		for (i in 1:9) { # 1-4 
			rect(4.3,1+(i-1)*(1/3),4.4,1+i*(1/3),col=colos[i],border=NA,xpd=NA)
			
		}
		for (i in 1:10) {
			text(4.8,1+(i-1)*(1/3),breakss[i],cex=2,xpd=NA)
		}
		text(4.25,4.3,expression(Deficit~(mm~d^-1)),cex=2,xpd=NA)
	}
}

dev.off()


