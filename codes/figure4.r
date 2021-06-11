source('default.R')
vars_name=c("evap","transpir","evapnu","inter")
te="newKsap5startDD11NOV0606gpsi1kf20copygmax1400day15"
evap_data=rep(NA,96*2*4)
dim(evap_data)=c(96,2,4)
for (ex in 1:2) {
	new=nc_open(paste0("CAX",expe[ex],te,"sechiba.nc"))
	for (vars in 1:4) {
		vardata=ncvar_get(new,vars_name[vars])
		mon=myfunction(vardata)
		evap_data[,ex,vars]=mon
	}
}


tiff("figure 4 evap partition 11June.tiff",width=900,height=540)

#dev.new(width=10,height=6)
colos=brewer.pal(5,"YlOrRd")
colos=c("white",colos)
colos=colos[6:1]
colos=brewer.pal(9,"YlOrRd")
colos=c(colos[9:1],'#FFFFFF')
par(cex=0.7, mai=c(0.1,0.05,0.2,0.1))
annual_prec=ann_prec
month_prec=myfunction(annual_prec)

for (ex in 1:2) { 
	new=nc_open(paste0("CAX",expe[ex],"newKsap5startDD11NOV0606gpsi1kf20copygmax1400day15sechiba.nc"))
	evap=ncvar_get(new,"evap")
	month_evap=myfunction(evap)
	if (ex==2) {
		month_prec[13:96]=month_prec[13:96]/2
	}
	deficit=month_prec-month_evap
	dim(deficit)=c(96,1)
	df=deficit
	df[df>0]=0.5
	df[df<(-3)]=-2.9 # 337
	par(fig=c((ex-1)*0.5+0.05,0.5+(ex-1)*0.5-0.012,0.68,0.75),new=TRUE)
	#image.plot(df,las=1,zlim=c(min(df),0),col=colos,cex.axis=2,xaxt="n",ylab="",xlab="",yaxt="n",axes=FALSE)
	image(df,las=1,zlim=c(-3,1),
	col=colos,cex.axis=2,xaxt="n",ylab="P-ET",xlab="",xaxs='i',yaxs='i',
	yaxt="n",axes=FALSE,cex.lab=2,breaks=c(-3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0,1),
	lab.breaks=c(-3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0,' '))
	# -3,-2.5,-2,-1.5,-0.5,-0.3,-0.2,-0.1,-0.01,0
	# old -3,-2.5,-2,-1.5,-1,-0.5,0
	axis(2,at=seq(-1,by=0.5,length.out=5),labels=seq(-1,by=0.5,length.out=5),
	cex.axis=2,col.axis="white",
	col="white",col.lab="white")
	axis(1, at=seq(0,by=1/8,length.out=10), 
			labels=seq(2001,by=1,length.out=10),cex.axis=2,col.axis="white")
}

wetp=c()

for (ex in 1:2) {
	par(fig=c((ex-1)*0.5,0.5+(ex-1)*0.5,0.1,0.8),new=TRUE)
	par(mar=c(5,5,3,2))
	evap=evap_data[,ex,1]
	tran=evap_data[,ex,2]
	soil=evap_data[,ex,3]
	canopy=evap_data[,ex,4]
	if (ex==1) {
		ylabs=expression(T~or~E~or~CE~(mm~d^-1))
	} else {
		ylabs=" "
	}
	plot(evap,type="l",col="#cbd5e8",ylim=c(0,6),xaxt="n",xlab="",
	main=expe[ex],cex.main=2,ylab=ylabs,
	cex.lab=2,cex.axis=2,yaxs="i",xaxs="i")
	##
	axis(2, at=seq(-1,7,by=1), labels = FALSE,cex.axis=2, lwd = 2)
	#axis(1,at=seq(0,6,by=1), labels = FALSE,lwd = 2)
	axis(1, at=seq(1,by=12,length.out=9), 
			labels=seq(2001,by=1,length.out=9),cex.axis=2,lwd=2)	
	axis(3,at=seq(1,by=12,length.out=9),cex.axis=2,lwd=2,tck=0,labels=FALSE)
	axis(4,at=seq(-1,7,by=1),tck=0,lwd=2,labels=FALSE)	
	##
		
	lines(soil,col="#b3e2cd")
	lines(soil+canopy,col="#fdcdac")
	xx=1:96
	polygon(c(xx,rev(xx)),c(rep(0,96),rev(soil)),col="#fdc086",border=NA)
	polygon(c(xx,rev(xx)),c(soil,rev(soil+canopy)),col="#7570b3",border=NA)
	polygon(c(xx,rev(xx)),c(soil+canopy,rev(evap)),col="#7fc97f",border=NA)
	
	#### emilieTEST0427MR001minus17keepDD0430less
	new=nc_open(paste0("CAX",expe[ex],
	"newKsap5startDD11NOV0606gpsi1kf20copygmax1400day15sechiba.nc"))
	evap=ncvar_get(new,"evap")
	month_evap=myfunction(evap)
	annual_prec=ann_prec
	month_prec=myfunction(annual_prec)
	if (ex==2) {
		month_prec[13:96]=month_prec[13:96]/2
	}
	deficit=month_prec-month_evap
	dim(deficit)=c(96,1)
	df=deficit
	df[df>0]=0
	
	if (ex==1) {
		for (year in 1:8) {
			xxdf=df[((year-1)*12+1):(year*12)]
			rows=which(xxdf<0)+(year-1)*12
			if (year<4) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==4 | year==5) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			
			if (year==7) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==8) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==6) {
				lines(c(rows[2]-0.5,rows[2]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
		}
	}
	if (ex==2) {
		for (year in 1:8) {
			xxdf=df[((year-1)*12+1):(year*12)]
			rows=which(xxdf<0)+(year-1)*12
			if (year<3) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==4 | year==7) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==3) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==5) {
				lines(c(rows[1]-0.5,rows[1]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==8) {
				lines(c(rows[2]-0.5,rows[2]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)]+0.5,rows[length(rows)]+0.5),c(-1,7),lty=2)
			}
			if (year==6) {
				lines(c(rows[3]-0.5,rows[3]-0.5),c(-1,7),lty=2)
				lines(c(rows[length(rows)-1]+0.5,rows[length(rows)-1]+0.5),c(-1,7),lty=2)
			}
		}
	}
	
	#arrows(x0 = 13, y0 = -1, x1 = 13, y1 = -0.3,xpd=TRUE,lwd=2)
	Arrows(x0 = 13, y0 = -1, x1 = 13, y1 = -0.3,xpd=TRUE,lwd=2,arr.width=0.15)
}

dev.off()