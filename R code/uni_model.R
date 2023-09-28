####################################
### updated: 09/26/23
### by Mikyoung Jun
####################################

NEED TO BE MODIFIED FURTHER AND COMMENTING

#### Note: this file contains codes to fit M1-3, but minor modification would be needed 
#### (in places commented below) to fit other M1-x models 

set.seed(0213)

### call libraries and necessary data ###

library("chron")
library("fields")
library("doMC")
registerDoMC(cores=28)

ncluster=28

### read data ###

popu

popu=as.matrix(popu)


## calculate PP (maximum population) for each year

PP=NULL

for(i in 3:14){

PP=c(PP, max(log(popu[,i]+1)))

}

popu1=NULL

for(i in 3:14){

popu1=cbind(popu1, log(popu[S.index,i]+1)/PP[i-2])
}


load("pop.ta.new.RData")

time1=data1[,3]##[events.ta[,1]==yy] ## julian day


pop.ta1=NULL
for(i in 1:length(pop.ta)){


pop.ta1=c(pop.ta1, log(pop.ta[i]+1)/PP[(data[i,1]-2001)])    ##[events.ta[,1]==yy]/sum(popu[,14])
}


load("pop.ta.new.RData")

time1=data1[,3]##[events.ta[,1]==yy] ## julian day

t_u=sort(unique(time1)) ## julian day

T0=min(t_u)
T1=max(t_u)

##load("afghan-spec-all.RData")
##data=data1

##load("afghan.RData")

##load("afg-jittered-data.RData")



### space-time triggering function ###

## below is an example of triggering function for M1-3

g=function(s,s_i, dt, theta, w, sigma,pp1,pp2,sigma1,SS){   
 
s=matrix(s, ncol=2) ## multiple rows
s_i=matrix(s_i, ncol=2) ## multiple rows matrix

d=rdist(s,s_i) ## Euclidean distance used (but chordal distance can be used instead)

ss=sigma+sigma1*mean(c(pp1,pp2))
temp1=theta*(1/w * exp(-dt/w))*1/(2*pi*SS^2) *exp(-d^2/(2*ss^2))

temp=sum(temp1)
return(temp)
}

### prep for numerical approximation of the integral of intensity ###

n_t=800 ## time resolution 
n_s=10000 ## spatial resolution 

S=as.matrix(popu[,1:2])
S.index=sort(sample(1:dim(S)[1], n_s))
S1=S[S.index,]
S1=as.matrix(S1)

area=62.71 ## approximate area of the country

s_area=area/n_s

tt=sort(runif(n_t, T0, T1)) ## time points for approximation

t_int=1/n_t


### loglikelihood function ###

logLik = function(par){
  
  print(par)

theta1=exp(par[1])/(1+exp(par[1]))

  w1=exp(par[2])

  sigma1=exp(par[3])

  mu0=par[5]
  mu1=0
  
sigma11=exp(par[4])-sigma1

SS=sigma1+sigma11##*max(popu1) ## !

a=floor(length(t_u)/ncluster)

g_i=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

temp=data1[time1==t_u[j],1:2]

if(j >1){

for(k in 1:(j-1)){

if(t_u[j]-t_u[k]<300){

temp1=data1[time1==t_u[k],1:2]

pp1=pop.ta1[time1==t_u[k]]
pp2=pop.ta1[time1==t_u[j]]

temp2=g(temp1, temp, (t_u[j]-t_u[k])/(T1-T0),theta1,w1,sigma1,pp1,pp2,sigma11,SS)

temp3=rbind(temp3, c(k,j,temp2))
}


}
}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(t_u))

g_i_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

temp=data1[time1==t_u[aa[i]],1:2]

for(k in 1:(aa[i]-1)){

##if(t_u[aa[i]]-t_u[k]<300){

temp1=data1[time1==t_u[k],1:2]

pp1=pop.ta1[time1==t_u[k]]
pp2=pop.ta1[time1==t_u[aa[i]]]

temp2=g(temp1, temp, (t_u[aa[i]]-t_u[k])/(T1-T0),theta1,w1,sigma1,pp1,pp2,sigma11,SS)

temp3=rbind(temp3, c(k,aa[i],temp2))
##}
}

return(temp3)
}

g_i=rbind(g_i, g_i_r)

##}

a=floor(length(tt)/ncluster)

g_tt=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

for(k in 1:length(t_u)){

if(t_u[k]<tt[j]){

##if(tt[j]-t_u[k]<300){

temp1=data1[time1==t_u[k],1:2]

pp1=pop.ta1[time1==t_u[k]]

tempp=month.day.year(tt[j])$year
temp2=g(temp1,S1, (tt[j]-t_u[k])/(T1-T0),theta1,w1,sigma1, pp1  ,popu1[,tempp-2002+1],sigma11,SS)

temp3=rbind(temp3, c(k,j,temp2))
##}

}
}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(tt))

g_tt_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

for(k in 1:length(t_u)){

if(t_u[k]<tt[aa[i]]){

##if(tt[aa[i]]-t_u[k]<300){

temp1=data1[time1==t_u[k],1:2]

pp1=pop.ta1[time1==t_u[k]]

tempp=month.day.year(tt[aa[i]])$year

temp2=g(temp1,S1, (tt[aa[i]]-t_u[k])/(T1-T0),theta1,w1,sigma1,pp1,popu1[,tempp-2002+1],sigma11,SS)

temp3=rbind(temp3, c(k,aa[i],temp2))
##}

}
}

return(temp3)

}


g_tt=rbind(g_tt, g_tt_r)

logL1=0

for(i in 2:length(t_u)){ 

##ith time event

temp=g_i[g_i[,2]==i,]
temp=matrix(temp, ncol=3)

temp=sum(temp[,3]) ## sum of g's up to ith event time

pp=pop.ta1[time1==t_u[i]]
temp=log((temp)+(exp(rep(mu0, length(pp)))))

temp=sum((temp))

logL1=logL1+temp
}

## second part (integral over space and time)

logL2=0

for(i in 1:length(tt)){

temp=g_tt[g_tt[,2]==i,]

temp=matrix(temp, ncol=3)

temp=sum(temp[,3])

tempp=month.day.year(tt[i])$year
pp1=popu1[,tempp-2002+1]

temp=(temp)+sum(exp(rep(mu0,length(pp1))))

temp=sum(temp)

logL2=logL2+temp
}

logL2=-logL2*t_int*s_area

##cat("log part1", logL1, " log part2", logL2,"\n")

temp=logL1+logL2

cat(temp, "final loglikelihood", "\n")

return(-temp) ## return negative likelihood for minimization

}


### numerical optimization ###

ini=c(2,-2.5,-2.5,0)

fit=optim(ini, logLik,control=list(maxit=10000000))

fit=nlm(logLik, ini, stepmax=10, print.level=2,hessian=TRUE)


