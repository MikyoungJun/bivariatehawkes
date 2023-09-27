####################################
### updated: 09/26/23
### by Mikyoung Jun
####################################

NEED TO BE MODIFIED FURTHER AND COMMENTING

#### Note: this file contains codes to fit M2-6, but minor modification would be needed (in places commented below) to fit other M2-x models 
####################### 

set.seed(0213)

load("Nigeria.RData") ## only Nigeria data saved (processed in my mac R)

library("fields")

library("doMC")
registerDoMC(cores=28)

ncluster=28

##install.packages("insol")


### new: due to problem with processing time now (change in machine?)

load("Nigeria.RData") ## only Nigeria data saved (processed in my mac R)

load("nigeria-processed.RData")


g=function(s,s_i, dt, theta, w, sigma,gamma){

s=matrix(s, ncol=2) ## multiple rows
s_i=matrix(s_i, ncol=2) ## multiple rows matrix

## treat time as scalar
##d=rdist.earth(s, s_i,miles=F)*1000 ## unit in meters
d=rdist(s, s_i) ## Euclidean distance

temp=sum(theta *1/w * exp(-dt/(T1-T0)/w)  *1/(2*pi*sigma^2*(1+dt/(T1-T0)/w)^(gamma)) *exp(-d^2/(2*sigma^2*(1+dt/(T1-T0)/w )^gamma )  ))
return(temp)
}

##N=20 ## integral bound for g12

g12=function(s,s_i, dt, theta, w, sigma,m1,m2){ ## alpha is the weight (??)

s=matrix(s, ncol=2) ## multiple rows
s_i=matrix(s_i, ncol=2) ## multiple rows matrix

s_i[,1]=s_i[,1]-m1
s_i[,2]=s_i[,2]-m2

## treat time as scalar
##d=rdist.earth(s, s_i,miles=F)*1000 ## unit in meters
d=rdist(s, s_i) ## Euclidean distance


temp=sum(theta *1/w * exp(-dt/(T1-T0)/w) *1 /(2*pi*sigma^2) * exp(-(d)^2/(2*sigma^2)) ) ## increase with spatial distance d
return(temp)
}



T0=max(min(data.b[,3]),min(data.f[,3]))
T1=min(max(data.b[,3]),max(data.f[,3]))

data.b=data.b[data.b[,3]>=T0 & data.b[,3]<=T1,]
data.f=data.f[data.f[,3]>=T0 & data.f[,3]<=T1,]

t_u.b=sort(unique(data.b[,3]))
t_u.f=sort(unique(data.f[,3]))


n_t=500 ## time resolution (at least 3000 needed...?)
n_s=ncluster*200 ## space resolution

##tt=seq(T0, T1,, n_t) ## regular "grid" for integration over time
tt=sort(runif(n_t, T0, T1))

load("finecov.RData")

index=sort(sample(1:dim(cov)[1], n_s))

##S=cbind(runif(n_s, 3,15), runif(n_s,5, 14)) ## irregular "grid" for integration over space

##S=matrix(S, ncol=2)

S=cov[index,1:2]
S=as.matrix(S)

t_int=1/n_t


area=0.06*0.048125*dim(cov)[1]
s_area=area/n_s


logLik = function(par){
  
  print(par)


### new parametrization!!

angle=(2*exp(par[1])/(1+exp(par[1]))-1)*pi/2
lambda1=0.1+exp(par[2])/(1+exp(par[2]))*(1-0.1)
lambda2=exp(par[3])/(1+exp(par[3]))*lambda1

A=matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)),2,2, byrow=FALSE)
B=diag(c(lambda1, lambda2),2,2)

A1=A %*% B %*% solve(A)

theta1=A1[1,1]
theta2=A1[2,2]

alpha1=A1[1,2]+exp(par[4])/(1+exp(par[4]))*(1-lambda1)

alpha2=A1[2,1]

cat(c(theta1, theta2, alpha1, alpha2), "\n")

A=matrix(c(theta1, alpha1, alpha2, theta2), 2,2, byrow=TRUE)
B=eigen(A)$values
print(B)


  w1=exp(par[5])
w2=exp(par[6])
w12=exp(par[7])

  sigma1=exp(par[8])
sigma2=exp(par[9])
sigma12=exp(par[10])


gamma1=exp(par[13])/(1+exp(par[13]))
gamma2=exp(par[14])/(1+exp(par[14]))

m1=(par[11])
m2=par[12]




mu0=0
mu01=0


## save g for event times first (boko haram)

print("g for boko")

data=data.b
t_u=t_u.b

a=floor(length(t_u)/ncluster)

g_i=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

temp=data[data[,3]==t_u[j],1:2]

if(j >1){

for(k in 1:(j-1)){

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1, temp, t_u[j]-t_u[k],theta1,w1,sigma1,gamma1)

temp3=rbind(temp3, c(k,j,temp2))
}
}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(t_u))

if(length(t_u)/ncluster!=a){

g_i_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

temp=data[data[,3]==t_u[aa[i]],1:2]

for(k in 1:(aa[i]-1)){

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1, temp, t_u[aa[i]]-t_u[k],theta1,w1,sigma1,gamma1)

temp3=rbind(temp3, c(k,aa[i],temp2))
}

return(temp3)
}

g_i=rbind(g_i, g_i_r)
}

## save g for general time (for time integral)

a=floor(length(tt)/ncluster)

g_tt=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

for(k in 1:length(t_u)){

if(t_u[k]<tt[j]){

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1,S, tt[j]-t_u[k],theta1,w1,sigma1,gamma1)

temp3=rbind(temp3, c(k,j,temp2))
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

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1,S, tt[aa[i]]-t_u[k],theta1,w1,sigma1,gamma1)

temp3=rbind(temp3, c(k,aa[i],temp2))
}

}

return(temp3)
}


g_tt=rbind(g_tt, g_tt_r)

g_i.b=g_i
g_tt.b=g_tt


### now fulani

print("g for fulani")

data=data.f
t_u=t_u.f

a=floor(length(t_u)/ncluster)

g_i=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

temp=data[data[,3]==t_u[j],1:2]

if(j >1){

for(k in 1:(j-1)){

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1, temp, t_u[j]-t_u[k],theta2,w2,sigma2,gamma2)

temp3=rbind(temp3, c(k,j,temp2))
}
}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(t_u))


if(length(t_u)/ncluster!=a){

g_i_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

temp=data[data[,3]==t_u[aa[i]],1:2]

for(k in 1:(aa[i]-1)){

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1, temp, t_u[aa[i]]-t_u[k],theta2,w2,sigma2,gamma2)

temp3=rbind(temp3, c(k,aa[i],temp2))
}

return(temp3)
}


g_i=rbind(g_i, g_i_r)
}

## save g for general time (for time integral)

a=floor(length(tt)/ncluster)

g_tt=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

for(k in 1:length(t_u)){

if(t_u[k]<tt[j]){

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1,S, tt[j]-t_u[k],theta2,w2,sigma2,gamma2)

temp3=rbind(temp3, c(k,j,temp2))
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

temp1=data[data[,3]==t_u[k],1:2]

temp2=g(temp1,S, tt[aa[i]]-t_u[k],theta2,w2,sigma2,gamma2)

temp3=rbind(temp3, c(k,aa[i],temp2))
}

}

return(temp3)
}


g_tt=rbind(g_tt, g_tt_r)

g_i.f=g_i
g_tt.f=g_tt


print("g for bf")

### now cross (bf first)

data=data.b
t_u=t_u.b

t_u_1=t_u.f
data1=data.f

a=floor(length(t_u)/ncluster)

g_i=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

temp=data[data[,3]==t_u[j],1:2]

##if(j >1){

for(k in 1:length(t_u_1)){

if(t_u_1[k]<t_u[j]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(temp1, temp, t_u[j]-t_u_1[k],alpha1,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,j,temp2))
}
}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(t_u))

if(length(t_u)/ncluster!=a){

g_i_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

temp=data[data[,3]==t_u[aa[i]],1:2]

for(k in 1:length(t_u_1)){

if(t_u_1[k]<t_u[aa[i]]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(temp1, temp, t_u[aa[i]]-t_u_1[k],alpha1,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,aa[i],temp2))
}
}

return(temp3)
}


g_i=rbind(g_i, g_i_r)
}

## save g for general time (for time integral)

a=floor(length(tt)/ncluster)

g_tt=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

for(k in 1:length(t_u_1)){

if(t_u_1[k]<tt[j]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(S,temp1, tt[j]-t_u_1[k],alpha1,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,j,temp2))
}

}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(tt))


g_tt_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

for(k in 1:length(t_u_1)){

if(t_u_1[k]<tt[aa[i]]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(S,temp1, tt[aa[i]]-t_u_1[k],alpha1,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,aa[i],temp2))
}

}

return(temp3)
}


g_tt=rbind(g_tt, g_tt_r)

g_i.bf=g_i
g_tt.bf=g_tt


## now cross (fb)

print("g for fb")

data=data.f
t_u=t_u.f

t_u_1=t_u.b
data1=data.b

a=floor(length(t_u)/ncluster)

g_i=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

temp=data[data[,3]==t_u[j],1:2]

##if(j >1){

for(k in 1:length(t_u_1)){

if(t_u_1[k]<t_u[j]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(temp1, temp, t_u[j]-t_u_1[k],alpha2,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,j,temp2))
}
}
}
return(temp3)
}


aa=((ncluster)*a+1):(length(t_u))


if(length(t_u)/ncluster!=a){

g_i_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

temp=data[data[,3]==t_u[aa[i]],1:2]

for(k in 1:length(t_u_1)){

if(t_u_1[k]<t_u[aa[i]]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(temp1, temp, t_u[aa[i]]-t_u_1[k],alpha2,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,aa[i],temp2))
}
}

return(temp3)
}


g_i=rbind(g_i, g_i_r)

}

## save g for general time (for time integral)

a=floor(length(tt)/ncluster)

g_tt=foreach(i=1:ncluster, .combine="rbind") %dopar%{

temp3=NULL

for(j in (i-1)*a+(1:a)){

for(k in 1:length(t_u_1)){

if(t_u_1[k]<tt[j]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(S,temp1, tt[j]-t_u_1[k],alpha2,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,j,temp2))
}

}
}
return(temp3)
}

aa=((ncluster)*a+1):(length(tt))


g_tt_r=foreach(i=1:length(aa), .combine="rbind") %dopar%{

temp3=NULL

for(k in 1:length(t_u_1)){

if(t_u_1[k]<tt[aa[i]]){

temp1=data1[data1[,3]==t_u_1[k],1:2]

temp2=g12(S,temp1, tt[aa[i]]-t_u_1[k],alpha2,w12,sigma12,m1,m2)

temp3=rbind(temp3, c(k,aa[i],temp2))
}

}

return(temp3)
}


g_tt=rbind(g_tt, g_tt_r)


g_i.fb=g_i
g_tt.fb=g_tt


#### likelihood part for boko haram first

print("likelihood calculation")

## first part of the likelihood function (not integral)

logL1=0

t_u=t_u.b
data=data.b

for(i in 2:length(t_u)){ ## correct starting from 2??

##ith time event

temp=g_i.b[g_i.b[,2]==i,]
tempp=g_i.bf[g_i.bf[,2]==i,]

temp=matrix(temp, ncol=3)
tempp=matrix(tempp, ncol=3)

temp=sum(temp[,3]) ## sum of g's up to ith event time

tempp=sum(tempp[,3])


temp1=data[data[,3]==t_u[i],1:2]

temp1=matrix(temp1, ncol=2)

##print(dim(temp1))

temp=(temp+tempp)##+exp(mu0)##+exp(mu0+mu1*(t_u[i]-T0)/T1+mu2*((t_u[i]-T0)/T1)^2+mu3*temp1[,1]+mu4*temp1[,2]+mu5*temp1[,1]*temp1[,2]) ## lambda_i

temp=log(sum(temp)+10^(-20))

logL1=logL1+temp
}

print("first part done")


## second part (integral over space and time)

logL2=0

for(i in 1:length(tt)){

temp=g_tt.b[g_tt.b[,2]==i,]

temp=matrix(temp, ncol=3)

temp=sum(temp[,3])

tempp=g_tt.bf[g_tt.bf[,2]==i,]
tempp=matrix(tempp, ncol=3)

tempp=sum(tempp[,3])

temp=(temp+tempp)##+exp(mu0)##+exp(mu0+mu1*(tt[i]-T0)/T1+mu2*((tt[i]-T0)/T1)^2+mu3*S[,1]+mu4*S[,2]+mu5*S[,1]*S[,2]) ## lambda

temp=sum(temp)

logL2=logL2+temp
}


logL2=-logL2*t_int*s_area

cat("log part1", logL1, " log part2", logL2,"\n")

logL1.b=logL1
logL2.b=logL2


## now likleihood part for fulani


## first part of the likelihood function (not integral)

logL1=0

t_u=t_u.f
data=data.f

for(i in 2:length(t_u)){ ## correct starting from 2??

##ith time event

temp=g_i.f[g_i.f[,2]==i,]
tempp=g_i.fb[g_i.fb[,2]==i,]

temp=matrix(temp, ncol=3)
tempp=matrix(tempp, ncol=3)

temp=sum(temp[,3]) ## sum of g's up to ith event time
tempp=sum(tempp[,3])

temp1=data[data[,3]==t_u[i],1:2]

temp1=matrix(temp1, ncol=2)

##print(dim(temp1))


if(temp+tempp<=0){

temp=10^(-300)


}

if(temp+tempp>0){



temp=(temp+tempp)#+exp(mu01)##+exp(mu01+mu11*(t_u[i]-T0)/T1+mu21*((t_u[i]-T0)/T1)^2+mu31*temp1[,1]+mu41*temp1[,2]+mu51*temp1[,1]*temp1[,2]) ## lambda_i

}

temp=log(sum(temp)+10^(-20))

logL1=logL1+temp
}


print("first part done")


## second part (integral over space and time)

logL2=0

for(i in 1:length(tt)){

temp=g_tt.f[g_tt.f[,2]==i,]

temp=matrix(temp, ncol=3)

temp=sum(temp[,3])

tempp=g_tt.fb[g_tt.fb[,2]==i,]
tempp=matrix(tempp, ncol=3)

tempp=sum(tempp[,3])

temp=(temp+tempp)##+exp(mu01)##+exp(mu01+mu11*(tt[i]-T0)/T1+mu21*((tt[i]-T0)/T1)^2+mu31*S[,1]+mu41*S[,2]+mu51*S[,1]*S[,2]) ## lambda

temp=sum(temp)

logL2=logL2+temp
}


logL2=-logL2*t_int*s_area

cat("log part1", logL1, " log part2", logL2,"\n")

logL1.f=logL1
logL2.f=logL2



temp=logL1.b+logL2.b+logL1.f+logL2.f


cat("loglikelihood", temp, "\n")

return(-temp) ## return negative likelihood for minimization

}



##fit=optim(ini, logLik,control=list(maxit=10000000))



##print(fit)

ini=c(0.5127791 , 1.1774306 , 2.3361057, 5.1848053 , 0.9516783,  0.5296093, -2.3333569, -5.5684722 ,-3.1044084,  0.5135355 , 1.7365658,  9.4819854, -0.7316783,0)

##fit=optim(ini, logLik, control=list(maxit=100000000))
load("m2-6-revision-result.RData")

ini=fit$estimate

fit=nlm(logLik, ini, stepmax=20, print.level=2,hessian=TRUE)

save(fit, file="m2-6-revision-result-1.RData")

print(fit)





