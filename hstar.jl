
using Distributions, LinearAlgebra, Random,KernelEstimator
using DataFrames,CSV

### Epanechnikov kernel function
function epkernel(x)
    temp = 0.75 .* (1 .- x .* x)
    temp[x .<= -1] .= 0.
    temp[x .>= 1] .= 0.
    temp
  end

### Gaussian kernel
function gausskernel(x)
    temp = 1/sqrt(2*π) .* exp.(- 0.5 .* x .^ 2)
    temp
  end


#univariate local linear regression
function locallinear(x, y, xgrid, h)
        # x: sorted data covariates,  y: data response, xgrid: grid points in the range of x
        # h:bandwidth
        newx = reshape(repeat(x, length(xgrid)), length(x),:)'
        # a matrix of size ngrid*length(x)
        xgridm = reshape(repeat(xgrid, length(x)), :, length(x))
        # a matrix of size of ngrid*length(x)
        kweights = reshape(gausskernel((newx .- xgridm) ./ h), :, length(x))
        u = fill(1, length(x))'
        sn0 = u * kweights'
        sn2 = diag((newx .- xgridm) .^ 2 * kweights')
        sn1 = diag((newx .- xgridm) * kweights')
        Sn0 = reshape(repeat(sn0', length(x)),:, length(x))
        Sn2 = reshape(repeat(sn2, length(x)), :, length(x))
        Sn1 = reshape(repeat(sn1, length(x)),:, length(x))
        W = (Sn2 .- (newx .- xgridm) .* Sn1) .* kweights
        W1 = (Sn0 .* (newx .- xgridm) .- Sn1) .* kweights
        wsum = W * u'
        wsum1 = W1 * u'

        fitted0 = (W * y) ./ wsum
        # regression function
        fitted1 = (W1 * y) ./ wsum
        #derivative
        # smoothing matrix only when grid points=data points
        #Sm <- (diag(1/wsum, nrow = length(x), ncol = length(x)) %*% W)
        #when xgrid=x, Sm gives the smoothing matrix
        #traceS <-sum( diag(W)/wsum)

        return fitted0, fitted1, W,W1
end



### Hstar matrix  (Huang & Chen, 2008)
### Epanechnikov kernel function, normalizing & symmetric & sum to one
#x: data of 1-d explanatory variable, xgrid: grid points to calculate the integral, h: bandwidth
function hstar(x, xgrid, h)
    gridlength = length(xgrid)
    datalength = length(x)
    gridint = fill(xgrid[2] - xgrid[1],gridlength)
    localproj = zeros(datalength, datalength)
    bigK = zeros(datalength, gridlength)
    for i in 1:datalength
      bigK[i,:] = (gausskernel(((x[i] .- xgrid) ./ h)) ./ h)
    end
    adjvect = 1 ./ (bigK * gridint)
    bigKK = (Diagonal(adjvect) * bigK)
    for i in 1:gridlength
      #print(i)
      kweight = Diagonal(bigKK[:,i])
      bigX = hcat(fill(1, datalength), (x .- xgrid[i]))
      localH = bigX * inv(bigX' * kweight * bigX) * bigX' *  kweight
      localproj = localproj + (kweight * localH)
    end
    localproj = localproj .* gridint[1]
    defree = sum(Diagonal(localproj))
    return localproj,defree
  end

### count the number of data points in neighborhood of  plus-minus h around every grid point
function counts(xgrid,x,h)
  countdata=zeros(length(xgrid),1)
  for i in 1:length(xgrid)
    countdata[i]=sum(abs.(x .- xgrid[i]) .<= h )
  end
  countdata
end


#when n=100 set.seed(20160815)
#random seed  when n=400
using Random, Distributions

Random.seed!(2016081500)


### function to generate data for Example 1
### sim =number of simulations, n=sample size, rho=correlation between x1 and x2
function randomdata(sim, n, rho,theta=1,err=1)  ### dimensions

d=2
  ### data points
Ngrid1=401
Ngrid2=401
  data1=zeros(n,sim)
  data2=zeros(n,sim)
  yrep=yi= zeros(n,sim)
  gridint1=zeros(Ngrid1,sim)
  gridint2=zeros(Ngrid2,sim)

  ### bandwidth
  bwth1=zeros(1,sim)

  bwth2=zeros(1,sim)

  ### error term
if err == 1
  td=Normal()
  error=rand(td,n,sim) #matrix(rnorm(n*sim,0,1),nrow=n,ncol=sim,byrow=FALSE)
elseif err == 2
  td=TDist(5)
  error=rand(td,n,sim)
elseif err == 3
  td=Chi(5)
  error=rand(td,n,sim)
elseif err == 4
  td=Chi(10)
  error=rand(td,n,sim)
end

  ####mean and variance
  ####mean and variance matrix for multivariate normal########
  mu=zeros(d,1)
  vmat=[1 rho;rho 1]

  z=zeros(n,d)
  ### initial values for covariates
  x=zeros(n,d)


  ### generate data
  for i in 1:sim

    counts1=zeros(Ngrid1,1)
    counts2=zeros(Ngrid2,1)


    ### using the smallest bandwdith 0.25 to count the number of data points in neighborhood
    h1=0.1
    h2=0.1
    ### grid points to calculate Hstar integral
    intgrid1=zeros(Ngrid1,1)
    intgrid2=zeros(Ngrid2,1)


    ### in each neighboodhood, at least 3  data points; otherwise, re-generate the data

    while  all(counts1 .>= 3)+all(counts2 .>= 3) != 2
      ### generate Multivariate Normal
      tmd=MvNormal(vmat)

      z=rand(tmd,n)'

      x= map(x ->  atan(x) / pi, z)

      m1= map(m1 -> 0.5-6*m1^2+3*m1^3,x[:,1])
      m2= map(m2 -> theta*sin(pi*m2), x[:,2])
      yi = m1 .+ m2 .+ error[:,i]

      ### grid points to calculate Hstar integral from min_xi to max_xi
      intgrid1=range(minimum(x[:,1]),maximum(x[:,1]),length=Ngrid1)
      intgrid2=range(minimum(x[:,2]),maximum(x[:,2]),length=Ngrid2)

      counts1=counts(intgrid1,x[:,1],h1)
      counts2=counts(intgrid2,x[:,2],h2)
    end
    #end of  while

    #store the simulated data

    data1[:,i]=x[:,1]
    data2[:,i]=x[:,2]
    bwth1[i]=h1
    bwth2[i]=h2
    yrep[:,i] = yi
gridint1[:,i] =intgrid1
gridint2[:,i]=intgrid2
  end
  #end of loop i

  return (data1,data2, bwth1,bwth2,yrep, gridint1,gridint2)

end

## end of function random.data ###

#############################################################################################
######### Main program ############################


### function for simplified smooth backfitting

function ssb(xmat,y,intgridmat,bw)

    # counts

    nobs=length(y)
    ngrid=size(intgridmat)[1]
    nvar=size(xmat)[2]

    # hstar matrices
    hmat=zeros(nobs,nobs,nvar)
    dfhmat=zeros(nvar)
    for j in 1:nvar
      temp=hstar(xmat[:,j],intgridmat[:,j],bw[j])
      hmat[:,:,j]=temp[1]
      dfhmat[j]=temp[2]
    end

# parametric projections
pjmat=zeros(nobs,nobs,nvar)
for j in 1:nvar
  temp=hcat(ones(nobs,1), xmat[:,j])
  pjmat[:,:,j]=temp*inv(temp'*temp)*temp'
end

txmat=hcat(ones(nobs,1),xmat)
pmat=txmat*inv(txmat'*txmat)*txmat'
Imat=Matrix(I,nobs,nobs)

pcoef= txmat\y

resy=(Imat-pmat)*y
#
oldmstar=zeros(nobs,nvar)
for j in 1:nvar
  oldmstar[:,j]=hmat[:,:,j]*resy
end
#
newmstar=oldmstar

for i in 1:100
  for j in 1:nvar
    temp= newmstar[:,1:end .!=j]
    newmstar[:,j]=hmat[:,:,j]*(resy-sum(temp,dims=2))
  end
    println("Iteration ", i)
    convthrvec=zeros(nvar)
     for j in 1:nvar
      convthrvec[j]=sqrt(sum((newmstar[:,j] .- oldmstar[:,j]) .^2))/sqrt(sum(newmstar[:,j] .^ 2))
     end
     convthr=maximum(convthrvec)
  if convthr <= 1e-6 && i > 1
    conv=true
    break
  else
    oldmstar=newmstar
  end

end

# adding parametric back
ppvec=zeros(nobs,nvar+1)
for j in 1:(nvar+1)
  ppvec[:,j]=txmat[:,j] .* pcoef[j]
end

# function m
m=zeros(nobs,nvar+1)
m[:,1]=ppvec[:,1]
for j in 1:nvar
  m[:,j+1]=ppvec[:,j+1] .+ newmstar[:,j]
end

predicted=ppvec[:,1] .+ sum(m,dims=2)

# local linear
l=zeros(ngrid,nvar+1)
l[:,1] = fill(ppvec[1,1],ngrid)
for j in 1:nvar
  tempy=vec(y-sum(m[:,1:end .!=j+1],dims=2))
  l[:,j+1]=locallinear(xmat[:,j],tempy,intgridmat[:,j],bw[j])[1]
end
#
rdict=Dict("mu" =>m, "beta" => l, "parmet" => ppvec, "fit" => predicted, "Hstar"=>hmat, "Gjmat"=>pjmat, "Gmat"=>pmat, "Xmat"=>txmat)
# return
return rdict #m,l, ppvec, predicted, hmat,pjmat,pmat,txmat

end

######
## smooth backfitting estimators and optimal bandwidth
#using KernelEstimator
function optbwd(xmat,y,intgridmat,ibw=0.1)
    # necessary
      nobs=length(y)
      ngrid=size(intgridmat)[1]
      nvar=size(xmat)[2]
#

bwd=fill(ibw,nvar)

for j in 1:nvar
  m=ssb(xmat,y,intgridmat,bwd)["mu"]
  tempy= vec(y-sum(m[:,1:end .!=j+1],dims=2))
  bwd[j]=bwlocallinear(xmat[:,j],tempy)
#  ld[:,j+1]=locallinear(xmat[:,j],tempy,xmat[:,j],bwd[j])[1]
end
# final fit
fitobj =ssb(xmat,y,intgridmat,bwd)
get!(fitobj,"bwd", bwd)

return fitobj

end

##
##
function linex(α,β,z)
  d= β/α^2 *(exp(α*z)-(1+α*z))
return d
end


#### conditional bootstrap

function cboot(data, ncboot, bwur,rind=2)
# extracting simulation 1 for conditional bootstrap
n,s=size(data[1])
y=data[5][:,1]
xmat=hcat(data[1][:,1],data[2][:,1])
bw=[data[3][1],data[4][1]]
intgridmat=hcat(data[6][:,1],data[7][:,1])


# unrestr model
  uresm= ssb(xmat,y,intgridmat,bwur)
  fitur=uresm["fit"]
  #bwur=uresm["bwd"]
  # restricted model data
  xmat_res=xmat[:, 1:end .!=rind]
  intgridmat_res=intgridmat[:,1:end .!=rind]

  bw_res=bwur[1:end .!= rind]

  resm=ssb(xmat_res,y,intgridmat_res,bw_res)
  fitr=resm["fit"]

  rss1=sum((y .- fitur) .^ 2)
  rss0=sum((y .- fitr) .^ 2)

# chi-square statistics
  lambda= n/2 * (rss0-rss1)/rss1
   Q1= n* sum(linex.(1e-6,1,fitur-fitr))/rss1
   Q2= n* sum(linex.(0.2,1,fitur-fitr))/rss1
   Q3= n* sum(linex.(0.5,1,fitur-fitr))/rss1
   Q4= n* sum(linex.(1,1,fitur-fitr))/rss1
# F statistcs
   #Gmd=xmat_res*inv(xmat_res'*xmat_res)*xmat_res'
   #Gd= uresm["Gjmat"][:,:,rind]
   #Igmd=(I(nobs)-Gmd)*Gd
   #Pgmd=Igmd*inv(Igmd'*Igmd)*Igmd'
   #Ig=I(nobs)-uresm["Gmat"]
   #C= Pgmd+ Ig*uresm["Hstar"][:,:,rind]+uresm["Hstar"][:,:,rind]*Ig-Ig*uresm["Hstar"][:,:,rind]*uresm["Hstar"][:,:,rind]*Ig
   #D=Ig-reduce(+, [uresm["Hstar"][:,:,i]*Ig+Ig*uresm["Hstar"][:,:,i]- Ig*uresm["Hstar"][:,:,i]*uresm["Hstar"][:,:,i]*Ig for i in 1:nvar],dims=1)[1]
   #E=Pgmd+ Ig*uresm["Hstar"][:,:,rind]+uresm["Hstar"][:,:,rind]*Ig
   #Fglr=((y'*C*y)/(y'*D*y))*(tr(D)/tr(C))

#
  epsi_ures=y .- fitur

  cen_epsi_ures=epsi_ures .- mean(epsi_ures)

  # loop through number of bootstraps
lambdastar=zeros(ncboot)
Q1star= zeros(ncboot)
Q2star= zeros(ncboot)
Q3star= zeros(ncboot)
Q4star= zeros(ncboot)
#
for i in 1:ncboot
  ystar = resm["fit"] .+ sample(cen_epsi_ures,n)
  # bootstrapped
  uresmstar= ssb(xmat,ystar,intgridmat,bwur)
  resmstar=ssb(xmat_res,ystar,intgridmat_res,bw_res)


  rss1star=sum((ystar .- uresmstar["fit"]) .^ 2)
  rss0star=sum((ystar .- resmstar["fit"]) .^ 2)

  lambdastar[i]= n/2 * (rss0star-rss1star)/rss1star
#
Q1star[i]= n* sum(linex.(1e-6,1,uresmstar["fit"]-resmstar["fit"]))/rss1star
Q2star[i]= n* sum(linex.(0.2,1,uresmstar["fit"]-resmstar["fit"]))/rss1star
Q3star[i]= n* sum(linex.(0.5,1,uresmstar["fit"]-resmstar["fit"]))/rss1star
Q4star[i]= n* sum(linex.(1,1,uresmstar["fit"]-resmstar["fit"]))/rss1star
end

cbstat=convert(DataFrame,hcat(lambdastar,Q1star, Q2star, Q3star, Q4star))

# calculating power for all the simulated datasets
bplambda=zeros(s)
bpq1=zeros(s)
bpq2=zeros(s)
bpq3=zeros(s)
bpq4=zeros(s)
#
bplambda[1]=mean(lambda .<= lambdastar)
bpq1[1]=mean(Q1 .<= Q1star)
bpq2[1]=mean(Q2 .<= Q2star)
bpq3[1]=mean(Q3 .<= Q3star)
bpq4[1]=mean(Q4 .<= Q4star)
#

for j in 2:s
  yj=data[5][:,j]
  xmatj=hcat(data[1][:,j],data[2][:,j])

  #intgridmatj=hcat(data[6][:,1],data[7][:,1])
  uresmj= ssb(xmatj,yj,intgridmat,bwur)
  fiturj=uresmj["fit"]
  #bwur=uresm["bwd"]
  # restricted model data
  xmat_resj=xmatj[:, 1:end .!=rind]
  #intgridmat_resj=intgridmat[:,1:end .!=rind]

  #bw_res=bwur[1:end .!= rind]

  resmj=ssb(xmat_resj,yj,intgridmat_res,bw_res)
  fitrj=resmj["fit"]

  rss1j=sum((yj .- fiturj) .^ 2)
  rss0j=sum((yj .- fitrj) .^ 2)

  lambdaj= n/2 * (rss0j-rss1j)/rss1j
   Q1j= n* sum(linex.(1e-6,1,fiturj-fitrj))/rss1j
   Q2j= n* sum(linex.(0.2,1,fiturj-fitrj))/rss1j
   Q3j= n* sum(linex.(0.5,1,fiturj-fitrj))/rss1j
   Q4j= n* sum(linex.(1,1,fiturj-fitrj))/rss1j
# p-values
bplambda[j]=mean(lambdaj .<= lambdastar)
bpq1[j]=mean(Q1j .<= Q1star)
bpq2[j]=mean(Q2j .<= Q2star)
bpq3[j]=mean(Q3j .<= Q3star)
bpq4[j]=mean(Q4j .<= Q4star)
#
end

# Power calculation
Onealplambda=mean(bplambda .<= 0.01)
Fivealplambda=mean(bplambda .<= 0.05)
Tenalplambda=mean(bplambda .<= 0.1)

Onealpbpq1=mean(bpq1 .<= 0.01)
Fivealpbpq1=mean(bpq1 .<= 0.05)
Tenalpbpq1=mean(bpq1 .<= 0.1)

Onealpbpq2=mean(bpq2 .<= 0.01)
Fivealpbpq2=mean(bpq2 .<= 0.05)
Tenalpbpq2=mean(bpq2 .<= 0.1)

Onealpbpq3=mean(bpq3 .<= 0.01)
Fivealpbpq3=mean(bpq3 .<= 0.05)
Tenalpbpq3=mean(bpq3 .<= 0.1)

Onealpbpq4=mean(bpq4 .<= 0.01)
Fivealpbpq4=mean(bpq4 .<= 0.05)
Tenalpbpq4=mean(bpq4 .<= 0.1)

# return power values
One=[Onealplambda Onealpbpq1 Onealpbpq2 Onealpbpq3 Onealpbpq4]
Five=[Fivealplambda Fivealpbpq1 Fivealpbpq2 Fivealpbpq3 Fivealpbpq4]
Ten=[Tenalplambda Tenalpbpq1 Tenalpbpq2 Tenalpbpq3 Tenalpbpq4]

return One, Five, Ten
end


##### end of the program

#### Simulation setting from Fan et al 2005

### sim =number of simulations, n=sample size, rho=correlation between x1 and x2
function randomdata1(sim, n, rho,theta=1,err=1)  ### dimensions

d=2
  ### data points
Ngrid1=401
Ngrid2=401
  data1=zeros(n,sim)
  data2=zeros(n,sim)
  yrep=yi= zeros(n,sim)
  gridint1=zeros(Ngrid1,sim)
  gridint2=zeros(Ngrid2,sim)

  ### bandwidth
  bwth1=zeros(1,sim)

  bwth2=zeros(1,sim)

  ### error term
if err == 1
  td=Normal()
  error=rand(td,n,sim) #matrix(rnorm(n*sim,0,1),nrow=n,ncol=sim,byrow=FALSE)
elseif err == 2
  td=TDist(5)
  error=rand(td,n,sim)
elseif err == 3
  td=Chi(5)
  error=rand(td,n,sim)
elseif err == 4
  td=Chi(10)
  error=rand(td,n,sim)
end

  ####mean and variance
  ####mean and variance matrix for multivariate normal########
  #mu=zeros(d,1)
  #vmat=[1 rho;rho 1]

  #z=zeros(n,d)
  ### initial values for covariates
  x=zeros(n,d)


  ### generate data
  for i in 1:sim

    counts1=zeros(Ngrid1,1)
    counts2=zeros(Ngrid2,1)


    ### using the smallest bandwdith 0.25 to count the number of data points in neighborhood
    h1=0.1
    h2=0.1
    ### grid points to calculate Hstar integral
    intgrid1=zeros(Ngrid1,1)
    intgrid2=zeros(Ngrid2,1)


    ### in each neighboodhood, at least 3  data points; otherwise, re-generate the data

    while  all(counts1 .>= 3)+all(counts2 .>= 3) != 2
      ### generate Multivariate Normal
      #tmd=MvNormal(vmat)

      #z=rand(tmd,n)'

      u1=rand(n) .-0.5
      u2=rand(n) .-0.5

      # transformation for correlation
       x[:,1] = u1 .+ 0.5 .* u2
       x[:,2] = 0.5 .* u1 .+ u2
      #x= map(x -> 3 * atan(x) / pi, z)

      m1= map(m1 -> 0.5-6*m1^2+3*m1^3,x[:,1])
      m2= map(m2 -> theta*sin(pi*m2), x[:,2])
      yi = m1 .+ m2 .+ error[:,i]

      ### grid points to calculate Hstar integral from min_xi to max_xi
      intgrid1=range(minimum(x[:,1]),maximum(x[:,1]),length=Ngrid1)
      intgrid2=range(minimum(x[:,2]),maximum(x[:,2]),length=Ngrid2)

      counts1=counts(intgrid1,x[:,1],h1)
      counts2=counts(intgrid2,x[:,2],h2)
    end
    #end of  while

    #store the simulated data

    data1[:,i]=x[:,1]
    data2[:,i]=x[:,2]
    bwth1[i]=h1
    bwth2[i]=h2
    yrep[:,i] = yi
gridint1[:,i] =intgrid1
gridint2[:,i]=intgrid2
  end
  #end of loop i

  return (data1,data2, bwth1,bwth2,yrep, gridint1,gridint2)

end











##
#### test code
### Assume equally-spaced grid points  of size 401 on [-1,1] for either domemsion for estimation
nobs=[100 200]
ed=[1 2 3 4]
theta=[0 0.2 0.4 0.6 0.8 1.0]
ncboot=5000
nsim=2000
#
incr=1
onemat=zeros(length(nobs)*length(ed)*length(theta),5)
fivemat=zeros(length(nobs)*length(ed)*length(theta),5)
tenmat=zeros(length(nobs)*length(ed)*length(theta),5)

# iterating over all the parameters
for i in 1:length(nobs)
for j in 1:length(ed)
for k in 1:length(theta)

data=randomdata(nsim,nobs[i],0.6,theta[k],ed[i])
y=data[5][:,2]
xmat=hcat(data[1][:,2],data[2][:,2])
bw=[data[3][2],data[4][2]]
intgridmat=hcat(data[6][:,2],data[7][:,2])

#
bwur=optbwd(xmat,y,intgridmat)["bwd"]
res=cboot(data,ncboot,bwur,2)
#
onemat[incr,:]=res[1]
fivemat[incr,:]=res[2]
tenmat[incr,:]=res[3]
#
end
end
end
# write the results to csv
CSV.write("one.csv",convert(DataFrame,onemat))
CSV.write("five.csv",convert(DataFrame,fivemat))
CSV.write("ten.csv",convert(DataFrame,tenmat))














## unrestricted model data
nobs=200
data=randomdata(10,nobs,0.6)
y=data[5][:,1]
xmat=hcat(data[1][:,1],data[2][:,1])
bw=[data[3][1],data[4][1]]
intgridmat=hcat(data[6][:,1],data[7][:,1])


## restricted model data
rind=2
xmat_res=xmat[:, 1:end .!=rind]
bw_res=bw[1:end .!=rind]
intgridmat_res=intgridmat[:,1:end .!=rind]


# unrestr model
uresm= ssb(xmat,y,intgridmat,bw)
fitur=uresm[1]
resm=ssb(xmat_res,y,intgridmat_res,bw_res)

rss1=sum((y .- uresm[3]) .^ 2)
rss0=sum((y .- resm[3]) .^ 2)

lambda= nobs/2 * log(rss0/rss1)

epsi_ures=y .- uresm[3]

cen_epsi_ures=epsi_ures .- mean(epsi_ures)

# loop through number of bootstraps

ystar = resm[3] .+ sample(cen_epsi_ures,nobs)

# bootstrapped
uresmstar= ssb(xmat,ystar,intgridmat,bw)
resmstar=ssb(xmat_res,ystar,intgridmat_res,bw_res)


rss1star=sum((ystar .- uresmstar[3]) .^ 2)
rss0star=sum((ystar .- resmstar[3]) .^ 2)

lambdastar= nobs/2 * log(rss0star/rss1star)


















#################################################################
### simulation setting

#try rho=0.0, 0.3, 0.5, and 0.7
#when n=100 set.seed(20160815)
sim=1000; n=400; rho=0.0;
set.seed(2016081500)
rd=randomdata(sim,n,0.0)


rd.data1=rd$data1   #dimension n*sim
rd.data2=rd$data2
rd.h1=rd$bwth1
rd.h2=rd$bwth2
m1=2*sin(pi*rd.data1)
m2=2*sin(pi*rd.data2)


### calculate the sample correlations between covariates
COR=matrix(0,2,2)
for (i in 1:sim){
  COR=COR+cor(cbind(rd.data1[,i],rd.data2[,i]))
}
COR/sim
#average correlation among sim=1000

### response
#y=m1+m2+e   #dimension n*sim

### sample mean of response
my=matrix(0,nrow=n,ncol=sim,byrow=FALSE)
y<- rd$yrep
for (i in 1:sim){
  my[,i]=matrix(rep(mean(y[,i]),n),nrow=n,ncol=1,byrow=FALSE)
}
### centering response
y0=y-my

L=matrix(rep(1/n,n^2),nrow=n,ncol=n)
I=diag(rep(1,n))

###################################
######## assign the value of bandwidth here #######################
rd.h1=matrix(0.45, nrow=1, ncol=sim, byrow=FALSE)
rd.h2= matrix(0.45, nrow=1, ncol=sim, byrow=FALSE)



### initial values for estimates at data points, new are estimates after updating
m2star.old=matrix(0,nrow=n,ncol=sim,byrow=FALSE)
m1star.old=matrix(0,nrow=n,ncol=sim,byrow=FALSE)
m2star.new=matrix(0,nrow=n,ncol=sim,byrow=FALSE)
m1star.new=matrix(0,nrow=n,ncol=sim,byrow=FALSE)

### initial values for estimates at grid points
Gm2star=matrix(0,nrow=Ngrid1,ncol=sim,byrow=FALSE)
Gm1star=matrix(0,nrow=Ngrid2,,ncol=sim,byrow=FALSE)

# later check if sum of m1star=0   sum of m2star=0
sm1star.new=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
sm2star.new=matrix(0,nrow=1,ncol=sim,byrow=FALSE)

### estimates for y at data points
y.hat.BA=matrix(0,nrow=n,ncol=sim,byrow=FALSE)

### mse for estimates at data points
mse.BA1=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
mse.BA2=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
mse.BAt=matrix(0,nrow=1,ncol=sim,byrow=FALSE)

# how many interations
iterative.times.1=matrix(0,nrow=1,ncol=sim,byrow=FALSE)

#Revise MSE for estimates at data points in  [-1,1]
Rmse.BA1=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
Rmse.BA2=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
Rmse.BAt=matrix(0,nrow=1,ncol=sim,byrow=FALSE)

# to count the number of data points with both x1 and x2 in [-1,1]
indtall <-matrix(0,nrow=1,ncol=sim,byrow=FALSE)


### mse at grid points
Gmse.BA1=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
Gmse.BA2=matrix(0,nrow=1,ncol=sim,byrow=FALSE)
Gmse.BAt=matrix(0,nrow=1,ncol=sim,byrow=FALSE)


truem1grid <- 2*sin(pi*Xgrid1)
truem2grid <- 2*sin(pi*Xgrid2)

#####################################################################################################################
################# simulation procedure starts
for (i in 1:sim){
  ### Hstar
  intgrid1=seq(min(rd.data1[,i]),max(rd.data1[,i]),(max(rd.data1[,i])-min(rd.data1[,i]))/(Ngrid1))
  intgrid2=seq(min(rd.data2[,i]),max(rd.data2[,i]),(max(rd.data2[,i])-min(rd.data2[,i]))/(Ngrid2))

  H1star=proj.ep.symmetric(rd.data1[,i],intgrid1, rd.h1[i])$Hstar
  H2star=proj.ep.symmetric(rd.data2[,i],intgrid2, rd.h2[i])$Hstar

  ### initial values
  m1star0=(H1star-L)%*%y[,i]
  m2star0=(H2star-L)%*%y[,i]


  m2star.old[,i]=m2star0
  m1star.old[,i]=(H1star-L)%*%(y[,i]-m2star.old[,i])

  ## for loop not needed, just set j=1 then j=count number of iterations
  for (j in 1:1){
    dis1=1
    while ((dis1>0.00001) && (j<50)) {
      m2star.new[,i]=m2star0-(H2star-L)%*%(m1star.old[,i])
      m1star.new[,i]=m1star0-(H1star-L)%*%(m2star.new[,i])

      dis1=max(sqrt(sum((m2star.new[,i]-m2star.old[,i])^2))/sqrt(sum((m2star.new[,i])^2)),
               sqrt(sum((m1star.new[,i]-m1star.old[,i])^2))/sqrt(sum((m1star.new[,i])^2)))
      m2star.old[,i]=m2star.new[,i]
      m1star.old[,i]=m1star.new[,i]
      j=j+1
    }
    #### end of while
    iterative.times.1[1,i]=j
  }

  # to check if the sum of m1star or m2star =0
  sm1star.new[i]=sum(m1star.new[,i]);  sm2star.new[i]=sum(m2star.new[,i])

  y.hat.BA[,i]=my[,i]+m1star.new[,i]+m2star.new[,i]

  ### calculate mse  at all data points
  mse.BA1[1,i]=sum((m1[,i]-m1star.new[,i])^2)/n
  mse.BA2[1,i]=sum((m2[,i]-m2star.new[,i])^2)/n
  mse.BAt[1,i]=sum( (m1[,i]-m1star.new[,i]+ m2[,i]-m2star.new[,i] )^2)/n

  ### calculate mse for estimates at data points on [-1,1]
  ind1<- abs(rd.data1[,i])<1
  ind2<- abs(rd.data2[,i])<1
  indt<- ind1*ind2
  indtall[1,i]= sum(indt)
  Rmse.BA1[1,i]=sum((m1[indt,i]-m1star.new[indt,i])^2)/sum(indt)
  Rmse.BA2[1,i]=sum((m2[indt,i]-m2star.new[indt,i])^2)/sum(indt)
  Rmse.BAt[1,i]= sum((m1[indt,i]-m1star.new[indt,i]  +m2[indt,i]-m2star.new[indt,i])^2) /sum(indt)


  #### calculate mse for estimates at grid points on [-1,1]
  tempy1<- y0[,i]- m2star.new[,i]
  Gm1star[,i]<- llr16ep2(rd.data1[,i], tempy1, Xgrid1, rd.h1[i])$m0

  tempy2<- y0[,i]- m1star.new[,i]
  Gm2star[,i]<- llr16ep2(rd.data2[,i], tempy2, Xgrid2, rd.h2[i])$m0

 Gmse.BA1[1,i]=sum((truem1grid -Gm1star[,i])^2)/Ngrid1
 Gmse.BA2[1,i]=sum((truem2grid -Gm2star[,i])^2)/Ngrid2

#assume Ngrid1=Ngrid2
 Gmse.BAt[1,i]=  sum( (truem1grid -Gm1star[,i] + truem2grid -Gm2star[,i])^2) /Ngrid1

  cat(i,"/")
}
######### end of i=1 to sim

print(rd.h1[1,1])
print(rho)

COR/sim

mean(iterative.times.1);sd(as.numeric(iterative.times.1)); range(iterative.times.1);

### mse for estimates at all data points
mean(mse.BA1);sd(as.numeric(mse.BA1))
mean(mse.BA2);sd(as.numeric(mse.BA2))
mean(mse.BAt);sd(as.numeric(mse.BAt))

### mse for estimates at data points in [-1,1]

mean(indtall);sd(as.numeric(indtall)); range(indtall);
mean(Rmse.BA1);sd(as.numeric(Rmse.BA1))
mean(Rmse.BA2);sd(as.numeric(Rmse.BA2))
mean(Rmse.BAt);sd(as.numeric(Rmse.BAt))

### mse for estimates at grid points
mean(Gmse.BA1);sd(as.numeric(Gmse.BA1))
mean(Gmse.BA2);sd(as.numeric(Gmse.BA2))
mean(Gmse.BAt);sd(as.numeric(Gmse.BAt))

mse.BA=round(matrix(c(mean(mse.BA1),sd(as.numeric(mse.BA1)),
                      mean(mse.BA2),sd(as.numeric(mse.BA2)),
                      mean(mse.BAt),sd(as.numeric(mse.BAt))
),nrow=2,ncol=3,byrow=FALSE),6);mse.BA


Rmse.BA=round(matrix(c(mean(Rmse.BA1),sd(as.numeric(Rmse.BA1)),
                      mean(Rmse.BA2),sd(as.numeric(Rmse.BA2)),
                      mean(Rmse.BAt),sd(as.numeric(Rmse.BAt))
),nrow=2,ncol=3,byrow=FALSE),6);Rmse.BA


print('grid points')
Gmse.BA=round(matrix(c(mean(Gmse.BA1),sd(as.numeric(Gmse.BA1)),
                      mean(Gmse.BA2),sd(as.numeric(Gmse.BA2)),
                      mean(Gmse.BAt),sd(as.numeric(Gmse.BAt))
),nrow=2,ncol=3,byrow=FALSE),6);Gmse.BA
