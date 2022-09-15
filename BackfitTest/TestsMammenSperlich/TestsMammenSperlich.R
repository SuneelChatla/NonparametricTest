##
##   Stefan Sperlich    23.03.2020
##    TestsMammenSperlich.R    calling dll from Fortran
##
##    first a test for constant 1D separable nonpar. function, say S_n
##      (eg. constant vs varying coeff., insignificant add. component)
##    then a test for additive 2D interaction, say T_n 
##
##    trimming works differently for them, see below
##    bandwidths h, hd to be set by hand
##    both tests consist of two dll calls of which
##      the 1st estimates H0 and H1 model and test statistic    
##      the 2nd performs golden cut bootstrap and calculates p-value
##    All based on smoothed backfitting procedures
## 
rm(list = ls());   
#install.packages(c("wsbackfit","np","mgcv","plot3D")) # when running the 1st time
library(KernSmooth); library(mgcv); library(np); library(wsbackfit);
require(graphics); library(plot3D);  Sys.setenv(LANG = "en")
#setwd("C:\\SPERLICH\\FORTI\\BACKFIT")  # set your own path !

###################################################################
#########  for test S_n  ##########################################

dyn.load("sbvctest.dll");  is.loaded("sbvctest")
dyn.load("sbvctestb.dll"); is.loaded("sbvctestb")
getLoadedDLLs()   # check if all 4 dlls are loaded and TRUE

#######  headers of varying coeff test FORTRAN SUBROUTINES ######## 
###################################################################
# subroutine sbvctest(x,z,y,n,d,h,mf,ng,maxit,conv,m,hd,to,m0)
#  INPUT
#   x(n,d),z(n,d),y(n),h(d),mf(ng,d) else scalars, n,d,ng,m integers
#      where m is number of gridpoints cut on each side in statistic 
#  OUTPUT :
#   to    = test statistic
#   y     = (n) estimate of E[Y|X,Z] under H_0 (without constant)
#   maxit = number of iterations used
#   mf    = function estimates on grid  
#   m0    = estimated constant (under H_1)
###################################################################
# subroutine sbvctestb(x,z,y,n,d,h,ng,mb,mit,cc,m,hd,to,nb,pv,rd) 
## x(n,d),z(n,d),y(n),h(d),mb(n) else scalars, n,d,ng,nb,m integers
## mb is hat{y} under H0 (without constant), cc is convcrit and mit is maxit
#  OUTPUT :
#   pv =  p-value of test

############## test for varying coeff estimation: sbvcest ##########
# set Integers variables
n=as.integer(150)       # Number of observations
d=as.integer(3)         # Number of dimensions
ng=as.integer(40)       # Grid length (not too small, used for integration)
maxit=as.integer(30)    # Maximal number of iterations
m=as.integer(2)         # for trimming in test, see above
nb=as.integer(500)      # number of bootstrap samples
# set Matrix/Array/Vector
x=matrix(runif(n*d,0,2),n,d)     # Design matrix X for nonparametric components
mf= matrix(0,ng,d)               # estimates on grid under H1 (output)
z=matrix(rnorm(n*d),n,d) + 2.0   # Design matrix Z, linear multipliers
h=sqrt(diag(var(x)))/n^(1/5)*1.5 # Bandwidth vector for estimation
# set Double precision
conv=as.numeric(0.0001)   # Convergence criteria in backfitting
hd=as.numeric(sd(x[,1])/n^(1/5))   # bandwidth for test statistic
to=as.numeric(999)        # value of test statistic
m0=as.numeric(0)          # output-constant of model (c in notation below)
rd=as.numeric(1000.0)    # random seed for replicability of simulations
pv=as.numeric(1.0)        # initialize p-value 
# simulate a varying coefficient model Y = sum_j b_j(x_j) z_j + error
y<-(sin(pi*x[,2])+2)*z[,2]+(sin(pi*x[,3])+2)*z[,3]+rnorm(n) 
y<- y + (2+0.0*sin(pi*x[,1]))*z[,1]  # b_1(x_1)= 2 + a sin(pi x_1)   
                                     # a neq 0 gives alternative
xg<-mf; for(j in 1:d){  # equispaced grid 
   s <- ( max(x[,j])-min(x[,j]) )/(ng-1)
   for(i in 1:ng){ xg[i,j] <- min(x[,j]) + (i - 1)*s }      }   
# call first routine, assign to list erg
erg<-.Fortran("sbvctest",x=x,z=z,y=y,n=n,d=d,h=h,mf=mf,xg=xg
           ,ng=ng,maxit=maxit,conv=conv,m=m,hd=hd,to=to,m0=m0)
# plot to see nonparametric function estimates
 #windows(8,4); par(mfrow = c(1,3),mar=c(3,3,3,3)); 
 #for(j in 1:3){plot(xg[,j],sin(pi*xg[,j])+2); lines(xg[,j],erg$mf[,j])} 
# call bootstrap routine - provide test statistic and H0 model
to<-erg$to; mb<-erg$y;             # mit<-erg$maxit
pt<-.Fortran("sbvctestb",x=x,z=z,y=y,n=n,d=d,h=h,ng=ng,mb=mb
         ,mit=maxit,cc=conv,m=m,hd=hd,to=to,nb=nb,pv=pv,rd=rd)
pt$pv  # the bootstrap p-value

# unload and clean memory
rm(list = ls()); dyn.unload("sbvctest.dll"); dyn.unload("sbvctestb.dll");  



###################################################################
#########  for test T_n  ##########################################

dyn.load("sbiatest.dll");  is.loaded("sbiatest")
dyn.load("sbiatestb.dll"); is.loaded("sbiatestb")
getLoadedDLLs()   # check if all 4 dlls are loaded and TRUE

#######  headers of add. interaction test FORTRAN SUBROUTINES ###### 
####################################################################
# subroutine sbiatest(x,n,d,h,hd,ng,y,to,trim,ra,m3,m12)  
## x(n,d),y(n),h(d),m3(ng,d-2),m12(ng*ng) else scalars, n,d,ng integers 
#    attention: all covariates X rescaled to about [-ra/2;ra/2]
#    trim gives absolute boundary size at each side of range  
#  OUTPUT :
#   y  = regression estimates under H_0 without constant
#   to = test statistic value
#   m3 = function estimates (for x3,... without x1,x2)  
#   m12= function estimates (for x1,x2 under H_0 decomposition) 
###################################################################
# subroutine sbiatestb(x,y,yb,n,d,h,hd,pv,ng,nb,rd,to,trim,ra) 
## x(n,d),y(n),yb(n),h(d) else scalars, n,d,ng,nb integers 
#   where yb is hat{y} under H0 (without constant), rd is random seed 
#    Attention: all covariates X rescaled to about [-ra/2;ra/2] 
#  OUTPUT :
#   pv =  p-value of test

############### test for additive interaction estimation #########
# Integers variables
n=as.integer(200)         # Number of observations
d=as.integer(3)           # Number of dimensions
ng=as.integer(25)         # Grid length, not too small as used for integration
nb=as.integer(250)        # number of bootstrap samples
# Matrix/Array/Vector
x=matrix(rnorm(n*d),n,d)          # Design matrix X 
hs=as.vector(sqrt(diag(var(x))))  # Vector of standard deviations
h<-hs*(100/n)^(1/6)               # Bandwidth vector
m3 = matrix(0,ng,d-2)             # additive 1D function estimates on grid
m12 = as.vector(rep(0,ng*ng))     # add.interaction function on grid
# Double precision
hd=as.numeric(0.5/n^(1/6)) # bandwidth for test
hd<-hd*hs[1:2]             # scaled for X
to=as.numeric(999)         # value of original test statistic
ra=as.numeric(4.6)         # considered range of covariates, see above 
trim=as.numeric(1.96)      # trimming (limits), see above
rd=as.numeric(1000.)       # random seed for replicability of simulations
pv=as.numeric(1.0)         # initialize p-value 
# simulate an additive interaction model Y = c + b12(x_1,x_2) + sum_j>2 b_j(x_j) + error
y<-2*x[,1]+2*sin(2*x[,2])+x[,3]^2+rnorm(n) #+x[,1]*x[,2] # interaction
# call first routine, assign to list erg
erg<-.Fortran("sbiatest",x=x,n=n,d=d,h=h,hd=hd,ng=ng,y=y
                        ,to=to,trim=trim,ra=ra,m3=m3,m12=m12)
# plot to see nonparametric function estimates
   xg = as.vector(rep(0,ng)); s <- ra / (ng-1)
   for(i in 1:ng){ xg[i] <- (i - 1.0)*s  - 0.5*ra }     
windows(8,4); par(mfrow = c(1,2),mar=c(3,3,3,3)); 
 mg<-erg$m3[,1]; mg<-mg-mean(mg) ; plot(xg,mg,type="l"); 
 mg<-xg^2; mg<-mg-mean(mg); lines(xg,mg,type="p"); 
  xg<-rep(xg,ng); xg<-cbind(xg,sort(xg))
scatter3D(xg[,1],xg[,2],erg$m12)
# call bootstrap routine - provide original test statistic and H0 model
to<-erg$to; mb<-erg$y
pt<-.Fortran("sbiatestb",x=x,y=y,yb=mb,n=n,d=d,h=h,hd=hd,pv=pv
               ,ng=ng,nb=nb,rd=rd,to=to,trim=trim,ra=ra)
pt$pv  # the bootstrap p-value

# clear memory and unload
rm(list = ls()); dyn.unload("sbiatest.dll"); dyn.unload("sbiatestb.dll");  

############ for saving you may use ... then don't do rm(list=ls())
#write(pvalue,file="resultS22w-2fN.txt",append = TRUE, sep = " ")
#write.table(x,file="out",append=F,sep=" ",na="NA",dec=".",row.names=T,col.names=T)

############### END #########################################

