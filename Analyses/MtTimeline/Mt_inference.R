# test
library(lubridate)
library(magrittr)
library(RColorBrewer)

cols= brewer.pal(n = 5, name = "PuBu")[-(1)]

setwd("C:/Users/Michael Lachmann/Documents/CV19")

pop=c(England=67e6,Belgium=12e6,Portugal=10e6,Spain=47e6)

# Multipliers for positive cases to real cases
mult=1/c(England=0.024,Belgium=0.06,Portugal=0.09,Spain=0.06)


a=read.csv("european_country_rt-plus-cases.csv",header = T,as.is=T)
a$date %<>% ymd


layout(rbind(c(1:2,5),c(3:4,5)),widths = c(0.4,0.4,0.15))
m=par("mar")
par(mar=c(2.1,6.1,4.1,2.1))
for( country in c("Belgium", "England","Spain","Portugal")) {
  b=a[a$country==country,]
  
  o=order(b$date,na.last = T)
  b=b[o,]

    b$inf=cumsum(b$cases)
  population =  pop[country]
  
  R0=max(b$mean_rt,na.rm=T)
  
  Rt=b$median_rt
  
  
  
  multiplier = mult[country]

  St=1 - b$inf * multiplier / population
  i=!is.na(Rt)
  Rt=Rt[i]
  St=St[i]
  D=b$date[i]
  
  
  S.star=0.2
  Mt = exp(R0 * ( log(St)/log(S.star) - (1-log(Rt)/log(R0))))
  par(cex.axis=1.6)
  par(cex.lab=2)
  par(cex.main=2)
  plot(D,Mt,log="y",ylim=c(0.001,1),xlab="Date",ylab="NPI transmission multiplier",main=country ,type="n")
  S.stars=1-c(0.7,0.6,0.4,0.2)

  for( i in seq_along(S.stars) ){
    S.star=S.stars[i]
    Mt = exp(R0 * ( log(St)/log(S.star) - (1-log(Rt)/log(R0))))
    lines( D, Mt,col=cols[i],lwd=2)  
  }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(mar=c(0,0,0,0))
legend("right",legend=paste("HIT=",rev(round(1-S.stars,2))),col=cols[1:4],lty=1,lwd=3,horiz = F,inset = 0,box.lty = 0,cex=2)
par(mar=m)
