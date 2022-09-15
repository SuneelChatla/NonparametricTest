
using Distributions, LinearAlgebra, Random,KernelEstimator
using DataFrames,CSV,Statistics,RCall


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
### Gaussian kernel function, normalizing & symmetric & sum to one
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



  ### Smoother matrix for Backfitting  (Opsomer, 200)
### Gaussian kernel function, normalizing & symmetric & sum to one
#x: data of 1-d explanatory variable, xgrid: grid points to calculate the integral, h: bandwidth
function Smat(x, xgrid, h)
  gridlength = length(xgrid)
  datalength = length(x)
  gridint = fill(xgrid[2] - xgrid[1],gridlength)
  localproj = zeros(datalength, datalength)
  bigK = zeros(datalength, gridlength)
  for i in 1:datalength
    bigK[i,:] = (gausskernel(((x[i] .- xgrid) ./ h)) ./ h)
  end
  #adjvect = 1 ./ (bigK * gridint)
  #bigKK = (Diagonal(adjvect) * bigK)
  for i in 1:gridlength
    #print(i)
    kweight = Diagonal(bigK[:,i])
    bigX = hcat(fill(1, datalength), (x .- xgrid[i]))
    localproj[i,:] = [1 0] * inv(bigX' * kweight * bigX) * bigX' *  kweight
  end
  defree = sum(Diagonal(localproj))
return localproj,defree
end

### count the number of data points in neighborhood of  plus-minus h around every grid point

function mcounts(xgrid,x,h)
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
function randomdata(sim, n, rho,theta=1,err=1,beta=0)  ### dimensions

d=4
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
errort=zeros(n,sim)
error=zeros(n,sim)

  ### error term
if err == 1
  td= Normal()
  error=rand(td,n,sim) #matrix(rnorm(n*sim,0,1),nrow=n,ncol=sim,byrow=FALSE)

elseif err == 2
  td=TDist(5)
  errort=rand(td,n,sim)
  for j in 1:sim
    error[:,j]=(errort[:,j] .- mean(errort[:,j])) ./ std(errort[:,j])
  end
elseif err == 3
  td=Chi(5)
  errort=rand(td,n,sim)
  for j in 1:sim
    error[:,j]=(errort[:,j] .- mean(errort[:,j])) ./ std(errort[:,j])
  end
elseif err == 4
  td=Chi(10)
  errort=rand(td,n,sim)
  for j in 1:sim
    error[:,j]=(errort[:,j] .- mean(errort[:,j])) ./ std(errort[:,j])
  end
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
      bmult = (1 + beta*sqrt(var(0.5 .- 6 .* x[:,1] .^2 .+ 3 .* x[:,1] .^3)))

      m1= map(m1 -> bmult*(0.5-6*m1^2+3*m1^3),x[:,1])
      m2= map(m2 -> theta*sin(pi*m2), x[:,2])
      yi = m1 .+ m2 .+ error[:,i]

      ### grid points to calculate Hstar integral from min_xi to max_xi
      intgrid1=range(minimum(x[:,1]),maximum(x[:,1]),length=Ngrid1)
      intgrid2=range(minimum(x[:,2]),maximum(x[:,2]),length=Ngrid2)

      counts1=mcounts(intgrid1,x[:,1],h1)
      counts2=mcounts(intgrid2,x[:,2],h2)
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

function ssb(xmat,y,intgridmat,bw,nit=200)

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

txmat=hcat(ones(nobs,1),xmat )  #.- mean(xmat,dims=1))
pmat=txmat*inv(txmat'*txmat)*txmat'
Imat=Matrix(I,nobs,nobs)


#pcoef= txmat\y

resy=(Imat-pmat)*y
#
oldmstar=zeros(nobs,nvar)
for j in 1:nvar
  oldmstar[:,j]=hmat[:,:,j]*resy
end
#
newmstar=copy(oldmstar)
ino=0

for i in 1:nit
  for j in 1:nvar
    temp= newmstar[:,1:end .!=j]
    newmstar[:,j]=hmat[:,:,j]*(resy-sum(temp,dims=2))
  end
    #println("Iteration ", i)
    convthrvec=zeros(nvar)
     for j in 1:nvar
      convthrvec[j]=sqrt(sum((newmstar[:,j] .- oldmstar[:,j]) .^2))/sqrt(sum(newmstar[:,j] .^ 2))
     end
     convthr=maximum(convthrvec)
     ino=i
  if convthr <= 1e-6 && i > 1
    conv=true
    break
  else
    oldmstar=copy(newmstar)
  end

end

# adding parametric back
txcmat=xmat .- mean(xmat,dims=1)

#pcmat=txcmat*inv(txcmat'*txcmat)*txcmat'

pcoef= txcmat\(y-sum(newmstar,dims=2))

ppvec=zeros(nobs,nvar+1)
ppvec[:,1] .= mean(y)

for j in 2:(nvar+1)
  ppvec[:,j]=txcmat[:,j-1] .* pcoef[j-1]
end

# function m
m=zeros(nobs,nvar+1)
m[:,1]=ppvec[:,1]
for j in 2:(nvar+1)
  m[:,j]=ppvec[:,j] .+ newmstar[:,j-1]
end
#
predicted= sum(m,dims=2)

# local linear
l=zeros(ngrid,nvar)
#l[:,1] = fill(ppvec[1,1],ngrid)
for j in 2:(nvar+1)
  tempy=vec(y -sum(m[:,1:end .!=j],dims=2))
  l[:,j-1]=locallinear(xmat[:,j-1],tempy,intgridmat[:,j-1],bw[j-1])[1]
end
#
rdict=Dict("mu" =>m, "beta" => l, "parmet" => ppvec, "fit" => predicted, "Hstar"=>hmat, "Gjmat"=>pjmat, "Gmat"=>pmat, "Xmat"=>txmat, "mustar" => newmstar, "iter" => ino)
# return
return rdict #m,l, ppvec, predicted, hmat,pjmat,pmat,txmat

end

##############################
### function to perform classical backfitting fan & jiang 2005, Opsomer 2000.
##############################
function kcb(xmat,y,bw,nit=200)

  # counts

  nobs=length(y)
  #ngrid=size(intgridmat)[1]
  nvar=size(xmat)[2]

  # hstar matrices
  Imat=Matrix(I,nobs,nobs)
  cmat=Imat-fill(1,(nobs,nobs)) ./ nobs
  hmat=zeros(nobs,nobs,nvar)
  dfhmat=zeros(nvar)
  for j in 1:nvar
    temp=Smat(xmat[:,j],xmat[:,j],bw[j])
    hmat[:,:,j]=cmat*temp[1]
    dfhmat[j]=temp[2]
  end


resy= y.-mean(y)
#
oldmstar=zeros(nobs,nvar)
for j in 1:nvar
oldmstar[:,j]=hmat[:,:,j]*resy
end
#
newmstar=copy(oldmstar)

for i in 1:nit
for j in 1:nvar
  temp= newmstar[:,1:end .!=j]
  newmstar[:,j]=hmat[:,:,j]*(resy-sum(temp,dims=2))
end
  #println("Iteration ", i)
  convthrvec=zeros(nvar)
   for j in 1:nvar
    convthrvec[j]=sqrt(sum((newmstar[:,j] .- oldmstar[:,j]) .^2))/sqrt(sum(newmstar[:,j] .^ 2))
   end
   convthr=maximum(convthrvec)
   #println("conv", convthr)
if convthr <= 1e-6 && i > 1
  conv=true
  break
else
  oldmstar=copy(newmstar)
end
end
#
predicted= sum(newmstar,dims=2) .+ mean(y)
#
rdict=Dict("mu" =>newmstar,  "fit" => predicted, "Smat"=>hmat)
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
  tempm=ssb(xmat,y,intgridmat,bwd)
  m=tempm["mu"]
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
  ystar =  resm["fit"]  .+ sample(cen_epsi_ures,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
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

      counts1=mcounts(intgrid1,x[:,1],h1)
      counts2=mcounts(intgrid2,x[:,2],h2)
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


## BostonHousing data
using DelimitedFiles

function housing(test=0.0;
                 file=joinpath(tempdir(),"housing.data"),
                 url="https://archive.ics.uci.edu/ml/machine-learning-databases/housing/housing.data")
    if !isfile(file)
        isdir(dirname(file)) || mkpath(dirname(file))
        @info("Downloading $url to $file")
        download(url, file)
    end
    data = readdlm(file)'
    # @show size(data) # (14,506)
    x = data[1:13,:]
    y = data[14:14,:]
    #x = (x .- mean(x,dims=2)) ./ std(x,dims=2) # Data normalization
    if test == 0
        xtrn = xtst = x
        ytrn = ytst = y
    else
        r = randperm(size(x,2))          # trn/tst split
        n = round(Int, (1-test) * size(x,2))
        xtrn=x[:,r[1:n]]
        ytrn=y[:,r[1:n]]
        xtst=x[:,r[n+1:end]]
        ytst=y[:,r[n+1:end]]
    end
    return (xtrn, ytrn, xtst, ytst)
end

##########################################
### Example 1: revised with 4 additive components
### function to generate data for Example 1
### sim =number of simulations, n=sample size, rho=correlation between x1 and x2

function randomdatarev(sim, n, rho,theta=1,err=1,beta=0)  ### dimensions are hard-coded: later do more automatic

  d=4
    ### data points
  Ngrid1=401
  #Ngrid2=401
    data1=zeros(n,sim)
    data2=zeros(n,sim)
    data3=zeros(n, sim)
    data4=zeros(n,sim)
    #
    m1mat=zeros(n,sim)
    m2mat=zeros(n,sim)
    m3mat=zeros(n, sim)
    m4mat=zeros(n,sim)
    #
    yrep=zeros(n,sim)
    yi=zeros(n,1)
    m1c=zeros(n,1)
    m2c=zeros(n,1)
    m3c=zeros(n,1)
    m4c=zeros(n,1)
    #
    gridint1=zeros(Ngrid1,sim)
    gridint2=zeros(Ngrid1,sim)
    gridint3=zeros(Ngrid1,sim)
    gridint4=zeros(Ngrid1,sim)
    ### bandwidth
    bwth1=zeros(1,sim)
    bwth2=zeros(1,sim)
    bwth3=zeros(1,sim)
    bwth4=zeros(1,sim)
    #
  errort=zeros(n,sim)
  error=zeros(n,sim)
  
    ### error term
  if err == 1
    td= Normal()
    error=rand(td,n,sim) #matrix(rnorm(n*sim,0,1),nrow=n,ncol=sim,byrow=FALSE)
  
  elseif err == 2
    td=TDist(5)
    errort=rand(td,n,sim)
    for j in 1:sim
      error[:,j]=(errort[:,j] .- mean(errort[:,j])) ./ std(errort[:,j])
    end
  elseif err == 3
    td=Chi(5)
    errort=rand(td,n,sim)
    for j in 1:sim
      error[:,j]=(errort[:,j] .- mean(errort[:,j])) ./ std(errort[:,j])
    end
  elseif err == 4
    td=Chi(10)
    errort=rand(td,n,sim)
    for j in 1:sim
      error[:,j]=(errort[:,j] .- mean(errort[:,j])) ./ std(errort[:,j])
    end
  end
  
    ####mean and variance
    ####mean and variance matrix for multivariate normal########
    mu=zeros(d,1)
    vmat= (1-rho)*Matrix(I,d,d)+rho*fill(1,(d,d))
  
    z=zeros(n,d)
    ### initial values for covariates
    x=zeros(n,d)
  
  
    ### generate data
  for i in 1:sim
  
      counts1=zeros(Ngrid1,1)
      counts2=zeros(Ngrid1,1)
      counts3=zeros(Ngrid1,1)
      counts4=zeros(Ngrid1,1)
  
  
      ### using the smallest bandwdith 0.25 to count the number of data points in neighborhood
      h1=0.1
      h2=0.1
      h3=0.1
      h4=0.1
      ### grid points to calculate Hstar integral
      intgrid1=zeros(Ngrid1,1)
      intgrid2=zeros(Ngrid1,1)
      intgrid3=zeros(Ngrid1,1)
      intgrid4=zeros(Ngrid1,1)
  
      ### in each neighboodhood, at least 3  data points; otherwise, re-generate the data
  
      while  all(counts1 .>= 5)+all(counts2 .>= 5)+all(counts3 .>=5)+ all(counts4 .>=5) != 4
        ### generate Multivariate Normal
        tmd=MvNormal(vmat)
  
        z=rand(tmd,n)'
  
        x= map(x ->  2atan(x) / pi, z)
        bmult = (1 + beta*sqrt(var(0.5 .-  x[:,1] .^2 .+ 3 .* x[:,1] .^3)))
  
        m1= map(m1 -> bmult*(0.5-m1^2+3*m1^3),x[:,1])
        m2= map(m2 -> theta*sin(pi*m2), x[:,2])
        m3= map(m3 -> m3*(1-m3), x[:,3])
        m4= map(m4 -> exp(2m4-1), x[:,4])
        
        
        m1c = m1 .- mean(m1)
        m2c = m2 .- mean(m2)
        m3c = m3 .- mean(m3)
        m4c = m4 .- mean(m4)

        yi = m1c .+ m2c .+ m3c .+ m4c .+ error[:,i]
  
        ### grid points to calculate Hstar integral from min_xi to max_xi
        intgrid1=range(minimum(x[:,1]),maximum(x[:,1]),length=Ngrid1)
        intgrid2=range(minimum(x[:,2]),maximum(x[:,2]),length=Ngrid1)
        intgrid3=range(minimum(x[:,3]),maximum(x[:,3]),length=Ngrid1)
        intgrid4=range(minimum(x[:,4]),maximum(x[:,4]),length=Ngrid1)
  
        counts1=mcounts(intgrid1,x[:,1],h1)
        counts2=mcounts(intgrid2,x[:,2],h2)
        counts3=mcounts(intgrid3,x[:,3],h3)
        counts4=mcounts(intgrid4,x[:,2],h4) 

      end
      #end of  while
  
      #store the simulated data
  
      data1[:,i]=x[:,1]
      data2[:,i]=x[:,2]
      data3[:,i]=x[:,3]
      data4[:,i]=x[:,4]
      #
      bwth1[i]=h1
      bwth2[i]=h2
      bwth3[i]=h3
      bwth4[i]=h4
      #
      yrep[:,i] = yi
  gridint1[:,i] =intgrid1
  gridint2[:,i]=intgrid2
  gridint3[:,i] =intgrid3
  gridint4[:,i]=intgrid4
  #
  m1mat[:,i]=m1c
  m2mat[:,i]=m2c
  m3mat[:,i]=m3c
  m4mat[:,i]=m4c
  #
   
end
    #end of loop i
  
return (data1,data2,data3, data4, bwth1,bwth2,bwth3, bwth4,yrep, gridint1,gridint2,gridint3,gridint4,m1mat,m2mat,m3mat,m4mat)
  
end
  
  ## end of function random.data.rev ###
  
 ###################################### 
  #### conditional bootstrap including Fan & Jiang's test with 4 covariates
####################################
  function cbootrev(data, ncboot, bwur,rind=2)
  # extracting simulation 1 for conditional bootstrap
  #rind=2
  n,s=size(data[1])
  y=data[9][:,1]
  xmat=hcat(data[1][:,1],data[2][:,1],data[3][:,1],data[4][:,1])
  bw=[data[5][1],data[6][1],data[7][1],data[8][1]]
  intgridmat=hcat(data[10][:,1],data[11][:,1],data[12][:,1],data[13][:,1])
  
  
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
  # GLR of Fan and jiang

  uresm_fj= kcb(xmat,y,bwur)
  fitur_fj=uresm_fj["fit"]
      
  resm_fj=kcb(xmat_res,y,bw_res)
  fitr_fj=resm_fj["fit"]
  
    rss1_fj=sum((y .- fitur_fj) .^ 2)
    rss0_fj=sum((y .- fitr_fj) .^ 2)
  
  # chi-square statistics of fan and jiang
    lambda_fj= n/2 * (rss0_fj-rss1_fj)/rss1_fj


  
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
  #
  epsi_ures_fj=y .- fitur_fj
  cen_epsi_ures_fj=epsi_ures_fj .- mean(epsi_ures_fj)
    # loop through number of bootstraps
  lambdastar=zeros(ncboot)
  Q1star= zeros(ncboot)
  Q2star= zeros(ncboot)
  Q3star= zeros(ncboot)
  Q4star= zeros(ncboot)
  lambdastar_fj=zeros(ncboot)
  #
for i in 1:ncboot
  #i=1
  Random.seed!(i)
    ystar =  resm["fit"]  .+ sample(cen_epsi_ures,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
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
# for fj 
ystar_fj =  resm_fj["fit"]  .+ sample(cen_epsi_ures_fj,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
# bootstrapped
uresmstar_fj= kcb(xmat,ystar_fj,bwur)
resmstar_fj=kcb(xmat_res,ystar_fj,bw_res)


rss1star_fj=sum((ystar_fj .- uresmstar_fj["fit"]) .^ 2)
rss0star_fj=sum((ystar_fj .- resmstar_fj["fit"]) .^ 2)

lambdastar_fj[i]= n/2 * (rss0star_fj-rss1star_fj)/rss1star_fj

end
  
  #cbstat=convert(DataFrame,hcat(lambdastar,Q1star, Q2star, Q3star, Q4star,lambdastar_fj))
  
  # calculating power for all the simulated datasets
  bplambda=zeros(s)
  bpq1=zeros(s)
  bpq2=zeros(s)
  bpq3=zeros(s)
  bpq4=zeros(s)
  bplambda_fj=zeros(s)
  #
  bplambda[1]=mean(lambda .<= lambdastar)
  bpq1[1]=mean(Q1 .<= Q1star)
  bpq2[1]=mean(Q2 .<= Q2star)
  bpq3[1]=mean(Q3 .<= Q3star)
  bpq4[1]=mean(Q4 .<= Q4star)
  bplambda_fj[1]=mean(lambda_fj .<= lambdastar_fj)
  #
  
  for j in 2:s
  #j=2 
  yj=data[9][:,j]
    xmatj=hcat(data[1][:,j],data[2][:,j], data[3][:,j],data[4][:,j])
  
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
  # for fj
  uresmstarj_fj= kcb(xmatj,yj,bwur)
resmstarj_fj=kcb(xmat_resj,yj,bw_res)


rss1starj_fj=sum((yj .- uresmstarj_fj["fit"]) .^ 2)
rss0starj_fj=sum((yj .- resmstarj_fj["fit"]) .^ 2)

lambdaj_fj= n/2 * (rss0starj_fj-rss1starj_fj)/rss1starj_fj
  
     # p-values
  bplambda[j]=mean(lambdaj .<= lambdastar)
  bpq1[j]=mean(Q1j .<= Q1star)
  bpq2[j]=mean(Q2j .<= Q2star)
  bpq3[j]=mean(Q3j .<= Q3star)
  bpq4[j]=mean(Q4j .<= Q4star)
  bplambda_fj[j]=mean(lambdaj_fj .<= lambdastar_fj)
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

  Onealplambda_fj=mean(bplambda_fj .<= 0.01)
  Fivealplambda_fj=mean(bplambda_fj .<= 0.05)
  Tenalplambda_fj=mean(bplambda_fj .<= 0.1)
  
  # return power values
  One=[Onealplambda Onealpbpq1 Onealpbpq2 Onealpbpq3 Onealpbpq4 Onealplambda_fj]
  Five=[Fivealplambda Fivealpbpq1 Fivealpbpq2 Fivealpbpq3 Fivealpbpq4 Fivealplambda_fj]
  Ten=[Tenalplambda Tenalpbpq1 Tenalpbpq2 Tenalpbpq3 Tenalpbpq4 Tenalplambda_fj]
  
  return One, Five, Ten
end
  
  
  ##### end of the program
  #####################################################################################
####### Function to compute the Backfitting test from Mammeen & Sperlich 2021 #######
####################################################################################

function SnM2021(xmat,y,bwur,rind)
  h=bwur
      # creatiin of dataframe
  df=DataFrame(xmat,["x$i" for i in 1:size(xmat)[2]])
  df[!,:y]=vec(y)
  #  using wsbackfit package to do the calcualtions
  R"library(wsbackfit)"
  R"m0=sback(formula = y ~ sb(x1, h = $h[1]) + sb(x2, h = $h[2]) + sb(x3, h=$h[3])+ sb(x4, h=$h[4]), $df)"
  R"effects=m0$effects"
  R"lincoef=m0$coeff"
  
  # return to julia
  @rget(effects)
  @rget(lincoef)
  # computation of combined effect (parametric + nonparametric)
  m=effects + xmat .* lincoef[2:length(lincoef)]'

  # density weights calcuation
  p2hat=kerneldensity(xmat[:,rind],h=h[rind])
  # sorting
  tmat=hcat(xmat[:,rind],m[:,rind],p2hat)
  tomat=tmat[sortperm(tmat[:,1]),:]
  # Riemann weights
  diffvec=zeros(length(xmat[:,rind]))
  diffvec[2:size(xmat)[1]]=[abs(tomat[i-1,1]-tomat[i,1])  for i in 2:size(xmat)[1]]
  # computing the statistic
  sn=sum((tomat[:,2] .^2) .* tomat[:,3] .* diffvec) 
  # return test statistic
  return sn, effects, lincoef
end

################################################################################
### Calculation of F-Statistics from models under H0 and H1 
###############################################################################
function Ftestcalc(resmodel, uresmodel,newy,rind)
  y=newy
  n=length(y)
  xmat=uresmodel["Xmat"][:,1:end .!= 1]
  Gmd= resmodel["Gmat"]  
  Gd= uresmodel["Gjmat"][:,:,rind]
  Igmd=(I(n)-Gmd)*xmat[:,rind] #Gd
  Pgmd=Igmd*inv(Igmd'*Igmd)*Igmd'
  Ig=I(n)-uresmodel["Gmat"]
  C= Pgmd+ Ig*uresmodel["Hstar"][:,:,rind]+uresmodel["Hstar"][:,:,rind]*Ig-Ig*uresmodel["Hstar"][:,:,rind]*uresmodel["Hstar"][:,:,rind]'*Ig #Pgmd
  D=Ig-reduce(+, [uresmodel["Hstar"][:,:,i]*Ig for i in 1:size(xmat)[2]],dims=1)[1]
  E=Pgmd + uresmodel["Hstar"][:,:,rind]*Ig
  # F test for GLR
  Fglr=((y'*C*y)/(y'*D*y))*(tr(D)/tr(C))
  # F test for LF
  EtE=E'*E
  Flft=((y'*EtE*y)/(y'*D*y))*(tr(D)/tr(EtE))
  # normalization constants
  #sn=0
  #dn=0
  #for i in 1:(n-1)
  #for j in (i+1):n
  #sn=sn+ C[i,j]^2
  #dn=dn+EtE[i,j]^2
  #end
  #end

  # normalized lambda
  #rk=tr(C)/sn
  # normalized q
  #sk=2*tr(EtE)/(2*dn)
  ## return
  return Fglr, Flft
end

######################################################################################
######## Conditional bootstrap: 
############    1. SSB GLR 
############    2. SSB LF 
############    3. CB GLR Fan & Jiang 2015
############    4. Backfitting Test  Mammen and Sperlich 2021
############    5. F tests
######################################################################################

function cbootrev_all(data, ncboot, bwur,rind=2)
  # extracting simulation 1 for conditional bootstrap
  #rind=2
  n,s=size(data[1])
  y=data[9][:,1]
  xmat=hcat(data[1][:,1],data[2][:,1],data[3][:,1],data[4][:,1])
  bw=[data[5][1],data[6][1],data[7][1],data[8][1]]
  intgridmat=hcat(data[10][:,1],data[11][:,1],data[12][:,1],data[13][:,1])
  
  ###################################
  ## For the GLR and LF based test statistics
  # unrestr model
  uresm= ssb(xmat,y,intgridmat,bwur)
  fitur=uresm["fit"]
    #bwur=uresm["bwd"]
    # restricted model data
  xmat_res=xmat[:, 1:end .!=rind]
  intgridmat_res=intgridmat[:,1:end .!=rind]
  bw_res=bwur[1:end .!= rind]
  # restricted model 
  resm=ssb(xmat_res,y,intgridmat_res,bw_res)
  fitr=resm["fit"]
  rss1=sum((y .- fitur) .^ 2)
  rss0=sum((y .- fitr) .^ 2)
  
  #  statistics
  lambda= n/2 * (rss0-rss1)/rss1
  Q1= n* sum(linex.(1e-6,1,fitur-fitr))/rss1
  Q2= n* sum(linex.(0.2,1,fitur-fitr))/rss1
  Q3= n* sum(linex.(0.5,1,fitur-fitr))/rss1
  Q4= n* sum(linex.(1,1,fitur-fitr))/rss1
  
  ########################################
  # GLR test of Fan and jiang
  uresm_fj= kcb(xmat,y,bwur)
  fitur_fj=uresm_fj["fit"]
      
  resm_fj=kcb(xmat_res,y,bw_res)
  fitr_fj=resm_fj["fit"]
  
  rss1_fj=sum((y .- fitur_fj) .^ 2)
  rss0_fj=sum((y .- fitr_fj) .^ 2)
  
  # chi-square statistics of fan and jiang
  lambda_fj= n/2 * (rss0_fj-rss1_fj)/rss1_fj

  #########################################
  ## F statistcs for GLR and LF tests
  #@infiltrate
  Fglr,Flft=Ftestcalc(resm,uresm,y,rind)
  ###############################################
  ## Backfit test based on Mammen & Sperlich 2021 Biometrica
  bfsn, npcoef, pcoef=SnM2021(xmat,y,bwur,rind)
  
  ###############################################
  ## Conditional bootstrap
  ###############################################
  # For GLR and LFT from Simplified
  epsi_ures=y .- fitur
  cen_epsi_ures=epsi_ures .- mean(epsi_ures)
  # For Fan and Jiang
  epsi_ures_fj=y .- fitur_fj
  cen_epsi_ures_fj=epsi_ures_fj .- mean(epsi_ures_fj)
  # For backfit test
  epsi_ures_ms= y .- pcoef[1] - xmat * pcoef[2:length(pcoef)] - sum(npcoef,dims=2)
  cen_epsi_ures_ms=  epsi_ures_ms .- mean(epsi_ures_ms) 
  # loop through number of bootstraps
  lambdastar=zeros(ncboot)
  Q1star= zeros(ncboot)
  Q2star= zeros(ncboot)
  Q3star= zeros(ncboot)
  Q4star= zeros(ncboot)
  lambdastar_fj=zeros(ncboot)
  #
  Fglrstar=zeros(ncboot)
  Flftstar=zeros(ncboot)
  #
  sn_ms=zeros(ncboot)
  #
  for i in 1:ncboot
  #i=1
      Random.seed!(i)
      ystar =  resm["fit"]  .+ sample(cen_epsi_ures,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
      # Bootstrapped GLR and LFT
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
      # F test
      Fstar1,Fstar2=Ftestcalc(resmstar,uresmstar,ystar,rind)       
      Fglrstar[i]=Fstar1[1]; Flftstar[i]=Fstar2[1]
      # for fj 
      ystar_fj =  resm_fj["fit"]  .+ sample(cen_epsi_ures_fj,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
      # bootstrapped
      uresmstar_fj= kcb(xmat,ystar_fj,bwur)
      resmstar_fj=kcb(xmat_res,ystar_fj,bw_res)
              
      rss1star_fj=sum((ystar_fj .- uresmstar_fj["fit"]) .^ 2)
      rss0star_fj=sum((ystar_fj .- resmstar_fj["fit"]) .^ 2)
      
      lambdastar_fj[i]= n/2 * (rss0star_fj-rss1star_fj)/rss1star_fj

      # For Backfit test
      ystar_ms=  pcoef[1] .+ xmat[:,1:end .!= rind] * pcoef[2:length(pcoef)][1:end .!= rind] + sum(npcoef[:, 1:end .!= rind],dims=2) .+ sample(cen_epsi_ures_ms,n) 
      sn_ms[i]=SnM2021(xmat,ystar_ms,bwur,rind)[1]
  end
  
  #cbstat=convert(DataFrame,hcat(lambdastar,Q1star, Q2star, Q3star, Q4star,lambdastar_fj))
  
  # calculating power for all the simulated datasets
  bplambda=zeros(s)
  bpq1=zeros(s)
  bpq2=zeros(s)
  bpq3=zeros(s)
  bpq4=zeros(s)
  #
  bpfglr=zeros(s)
  bpflft=zeros(s)
  #
  bplambda_fj=zeros(s)
  #
  bpms=zeros(s)
  #
  bplambda[1]=mean(lambda .<= lambdastar)
  bpq1[1]=mean(Q1 .<= Q1star)
  bpq2[1]=mean(Q2 .<= Q2star)
  bpq3[1]=mean(Q3 .<= Q3star)
  bpq4[1]=mean(Q4 .<= Q4star)
  #
  bpfglr[1]=mean(Fglr .<= Fglrstar)
  bpflft[1]=mean(Flft .<= Flftstar)
  #
  bplambda_fj[1]=mean(lambda_fj .<= lambdastar_fj)
  #
  bpms[1]=mean(bfsn .<= sn_ms)
  #
  for j in 2:s
      #j=2 
      yj=data[9][:,j]
      xmatj=hcat(data[1][:,j],data[2][:,j], data[3][:,j],data[4][:,j])
      
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
      #
          ## F statistcs for GLR and LF tests
      Fglrj,Flftj=Ftestcalc(resmj,uresmj,yj,rind)
      ###############################################
      ## Backfit test based on Mammen & Sperlich 2021 Biometrica
      bfsnj, npcoefj, pcoefj=SnM2021(xmatj,yj,bwur,rind)
      # for fj
      uresmstarj_fj= kcb(xmatj,yj,bwur)
      resmstarj_fj=kcb(xmat_resj,yj,bw_res)
      
      
      rss1starj_fj=sum((yj .- uresmstarj_fj["fit"]) .^ 2)
      rss0starj_fj=sum((yj .- resmstarj_fj["fit"]) .^ 2)
      
      lambdaj_fj= n/2 * (rss0starj_fj-rss1starj_fj)/rss1starj_fj
  
      # p-values
      bplambda[j]=mean(lambdaj .<= lambdastar)
      bpq1[j]=mean(Q1j .<= Q1star)
      bpq2[j]=mean(Q2j .<= Q2star)
      bpq3[j]=mean(Q3j .<= Q3star)
      bpq4[j]=mean(Q4j .<= Q4star)
      #
      bpfglr[j]=mean(Fglrj .<= Fglrstar)
      bpflft[j]=mean(Flftj .<= Flftstar)
      #
      bplambda_fj[j]=mean(lambdaj_fj .<= lambdastar_fj)
      #
      bpms[j]=mean(bfsnj .<= sn_ms)
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

  #
  Onealpfglr=mean(bpfglr .<= 0.01)
  Fivealpfglr=mean(bpfglr .<= 0.05)
  Tenalpfglr=mean(bpfglr .<= 0.1)

  Onealpflft=mean(bpflft .<= 0.01)
  Fivealpflft=mean(bpflft .<= 0.05)
  Tenalpflft=mean(bpflft .<= 0.1)
  
  Onealplambda_fj=mean(bplambda_fj .<= 0.01)
  Fivealplambda_fj=mean(bplambda_fj .<= 0.05)
  Tenalplambda_fj=mean(bplambda_fj .<= 0.1)
  
  Onealpms=mean(bpms .<= 0.01)
  Fivealpms=mean(bpms .<= 0.05)
  Tenalpms=mean(bpms .<= 0.1)
  # return power values
  One=[Onealplambda Onealpbpq1 Onealpbpq2 Onealpbpq3 Onealpbpq4 Onealplambda_fj Onealpfglr Onealpflft Onealpms]
  Five=[Fivealplambda Fivealpbpq1 Fivealpbpq2 Fivealpbpq3 Fivealpbpq4 Fivealplambda_fj Fivealpfglr Fivealpflft Fivealpms]
  Ten=[Tenalplambda Tenalpbpq1 Tenalpbpq2 Tenalpbpq3 Tenalpbpq4 Tenalplambda_fj Tenalpfglr Tenalpflft Tenalpms]
  
  return One, Five, Ten
end
  
  
  ##### end of the program
  




