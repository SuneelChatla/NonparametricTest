



##### H star for p=0
######################################################
function hstartest(x, xgrid, h, xtest,kindex,pindex)
  gridlength = length(xgrid)
  datalength = length(x)
  testlength=length(xtest)
  gridint = fill(xgrid[2] - xgrid[1],gridlength)
  localproj = zeros(datalength, datalength)
  localprojtest = zeros(testlength, datalength)
  bigK = zeros(datalength, gridlength)
  # test
  bigKtest=zeros(testlength, gridlength)
  #
  for i in 1:datalength
    if kindex==1
      bigK[i,:] = (gausskernel(((x[i] .- xgrid) ./ h)) ./ h)
    else
      bigK[i,:] = (epkernel(((x[i] .- xgrid) ./ h)) ./ h)
    end
  end
  for i in 1:testlength
    if kindex==1
      bigKtest[i,:] = (gausskernel(((xtest[i] .- xgrid) ./ h)) ./ h)
    else
      bigKtest[i,:] = (epkernel(((xtest[i] .- xgrid) ./ h)) ./ h)
    end
  end

  adjvect = 1 ./ (bigK * gridint)
  bigKK = (Diagonal(adjvect) * bigK)
  #
  adjvecttest = 1 ./ (bigKtest * gridint)
  bigKKtest = (Diagonal(adjvecttest) * bigKtest)
  #
  for i in 1:gridlength
    #print(i)
    kweight = Diagonal(bigKK[:,i])
    if pindex==0
      bigX = fill(1, datalength)
    elseif pindex==1
      bigX=hcat(fill(1, datalength), (x .- xgrid[i]))
    end
    localH = bigX * inv(bigX' * kweight * bigX) * bigX' *  kweight
    localproj = localproj + (kweight * localH)
    #
    kweighttest = Diagonal(bigKKtest[:,i])
    if pindex==0
      bigXtest = fill(1, testlength)
    elseif pindex==1
      bigXtest= hcat(fill(1, testlength), (xtest .- xgrid[i]))
    end
    localHtest = bigXtest * inv(bigX' * kweight * bigX) * bigX' *  kweight
    localprojtest = localprojtest + (kweighttest * localHtest)
  end
  localproj = localproj .* gridint[1]
  localprojtest = localprojtest .* gridint[1]
  
  #defree = sum(Diagonal(localproj))
  return localproj,localprojtest
end




include("ssb-inf.jl")

n=100
nsim=2
err=[1 2 3 4]

i=1

data=randomdatarev(nsim,n,0.6,1,err[i])
y=data[9][:,2]
xmat=hcat(data[1][:,2],data[2][:,2],data[3][:,2],data[4][:,2])
bw=[data[5][2],data[6][2],data[7][2],data[8][2]]
intgridmat=hcat(data[10][:,2],data[11][:,2],data[12][:,2],data[13][:,2])

#####################################################
######  Eigenvales for different bandwidths (1- sample)
#######################################################

x=vec(xmat[:,1])
y=vec(y[:,1])
xgrid=vec(intgridmat[:,1])
# h-grid
hrange=[0.1,1]
nhgrid=40
hgrid=range(hrange[1],hrange[2],length=nhgrid)

evnw=zeros(n,nhgrid)
evhst0=zeros(n,nhgrid)
evhst1=zeros(n,nhgrid)
P0mat=ones(n)*inv(ones(n)'*ones(n))*ones(n)'
P1mat=hcat(ones(n),x)*inv(hcat(ones(n),x)'*hcat(ones(n),x))*hcat(ones(n),x)'
evspl=zeros(n,nhgrid)
kindex=1

for i in 1:nhgrid
    evnw[:,i]=real(eigvals(hstartest(x, xgrid, hgrid[i], x,kindex,0)[1]-P0mat))
    evhst0[:,i]=real(eigvals(hstartest(x, xgrid, hgrid[i], x,kindex,0)[1]))
    evhst1[:,i]=real(eigvals(hstartest(x, xgrid, hgrid[i], x,kindex,1)[1]))
    evspl[:,i]=real(eigvals(hstartest(x, xgrid, hgrid[i], x,kindex,1)[1]-P1mat))
end

using Plots, LaTeXStrings


p1=plot(hgrid,evnw', legend=false,xlabel="Bandwidth(h)",title=L"H^*_0-P_1",xguidefontsize=14,yguidefontsize=14, linestyle=:auto);
p2=plot(hgrid,evhst0', legend=false,xlabel="Bandwidth(h)",ylabel="Eigenvalue",title=L"H^*_0",xguidefontsize=14,yguidefontsize=14, linestyle=:auto);
p3=plot(hgrid,evhst1', legend=false,xlabel="Bandwidth(h)",title=L"H^*_1",xguidefontsize=14,yguidefontsize=14, linestyle=:auto);
p4=plot(hgrid,evspl', legend=false,xlabel="Bandwidth(h)",title=L"H^*_1-P_{[1 ~ x]}",xguidefontsize=14,yguidefontsize=14, linestyle=:auto);

# Plots for Eigenvalues
plot(p2,p1,p3,p4, layout=(1,4),margin=5Plots.mm);
plot!(size=(1250,350));
savefig("eig-h-testing.pdf")











include("ssb-inf.jl")

##
#### test code
### Assume equally-spaced grid points  of size 401 on [-1,1] for either domemsion for estimation
nobs=[100 200]
ed=[1 2 3 4]
theta=[0 0.2 0.4 0.6 0.8 1.0]
ncboot=3000
nsim=1000
#
incr=1
onemat=zeros(length(nobs)*length(ed)*length(theta),5)
fivemat=zeros(length(nobs)*length(ed)*length(theta),5)
tenmat=zeros(length(nobs)*length(ed)*length(theta),5)

# iterating over all the parameters
for i in 1:length(nobs)
bwur=[0.10 0.10]
for j in 1:length(ed)
for k in 1:length(theta)

data=randomdata(nsim,nobs[i],0.6,theta[k],ed[i])
y=data[5][:,2]
xmat=hcat(data[1][:,2],data[2][:,2])
bw=[data[3][2],data[4][2]]
intgridmat=hcat(data[6][:,2],data[7][:,2])

#
#bwur= [0.15 0.15] #optbwd(xmat,y,intgridmat)["bwd"]
res=cboot(data,ncboot,bwur,2)
#
onemat[incr,:]=res[1]
fivemat[incr,:]=res[2]
tenmat[incr,:]=res[3]
global incr=incr+1
#
end
end
end
# write the results to csv
CSV.write("ones.csv",convert(DataFrame,onemat))
CSV.write("fives.csv",convert(DataFrame,fivemat))
CSV.write("tens.csv",convert(DataFrame,tenmat))

####
#latex export
using LatexPrint
thetavec=repeat(theta',length(nobs)*length(ed))
edvec=repeat(vcat(repeat(ed,length(theta))...),2)
samplevec=vcat(repeat(nobs,24)...)

#
myfivemat=zeros(48,7)
myfivemat=hcat(edvec,thetavec,fivemat)
tabular(myfivemat)

#
mytenmat=zeros(48,7)
mytenmat=hcat(edvec,thetavec,tenmat)
tabular(mytenmat)

####
##  Figure 1
### Wilks phenomenon
include("ssb-inf.jl")
#
n= 200
beta=[-1.5 0 1.5]
nsim=1000
rind=2
bwdcom=zeros(3,2)
#

b1glr=zeros(nsim,length(beta)*3)
rkglr=zeros(nsim,length(beta)*3)
b3glr=zeros(nsim,length(beta)*3)
#
b1lft1=zeros(nsim,length(beta)*3)
b1lft2=zeros(nsim,length(beta)*3)
b1lft3=zeros(nsim,length(beta)*3)
b1lft4=zeros(nsim,length(beta)*3)


sklft=zeros(nsim,length(beta)*3)


# iterating over all the parameters
incr=1
#
for i in 1:length(beta)

data=randomdata(nsim,n,0.6,1,1,beta[i])
y=data[5][:,2]
xmat=hcat(data[1][:,2],data[2][:,2])
bw=[data[3][2],data[4][2]]
intgridmat=hcat(data[6][:,2],data[7][:,2])
#
# unrestr model
  uresm=optbwd(xmat,y,intgridmat)
  bw=uresm["bwd"]
bwdcom[2,:]=bw'
bwdcom[1,:]=[bw[1]*1/3,bw[2]]
bwdcom[3,:]=[bw[1]*1.5,bw[2]]


#

for k in 1:3
#
bwur=bwdcom[k,:]
for j in 1:nsim

y=data[5][:,j]
xmat=hcat(data[1][:,j],data[2][:,j])
intgridmat=hcat(data[6][:,j],data[7][:,j])
#
# unrestr model
  uresm=ssb(xmat,y,intgridmat,bwur)
#  bwur=uresm["bwd"]
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


Gmd= resm["Gmat"]  #xmat_res*inv(xmat_res'*xmat_res)*xmat_res'
Gd= uresm["Gjmat"][:,:,rind]
Igmd=(I(n)-Gmd)* xmat[:,rind]
Pgmd=Igmd*inv(Igmd'*Igmd)*Igmd'
Ig=I(n)-uresm["Gmat"]
C= Pgmd+ Ig*uresm["Hstar"][:,:,rind]+uresm["Hstar"][:,:,rind]*Ig-Ig*uresm["Hstar"][:,:,rind]*uresm["Hstar"][:,:,rind]'*Ig #Pgmd
D=Ig-reduce(+, [uresm["Hstar"][:,:,i]*Ig+Ig*uresm["Hstar"][:,:,i]- Ig*uresm["Hstar"][:,:,i]*uresm["Hstar"][:,:,i]*Ig for i in 1:2],dims=1)[1]
E=Pgmd+ Ig*uresm["Hstar"][:,:,rind]+uresm["Hstar"][:,:,rind]*Ig
#Fglr=((y'*C*y)/(y'*D*y))*(tr(D)/tr(C  Pgmd
EtE=E'*E

sn=0
dn=0
for i in 1:(n-1)
for j in (i+1):n
sn=sn+ C[i,j]^2
dn=dn+EtE[i,j]^2
end
end

# normalized lambda
rk=tr(C)/(2*sn)
# normalized q
sk=tr(EtE)/(2*dn)

rkglr[j,incr]=rk
sklft[j,incr]=sk


# storing normalized lambda and q statistics
b1glr[j,incr]=lambda  #rk*
b1lft1[j,incr]=Q1   #sk*
b1lft2[j,incr]=Q2
b1lft3[j,incr]=Q3
b1lft4[j,incr]=Q4


# for j
end
# for k
incr=incr+1
end
# for i
end

##
using KernelEstimator
using   Plots

td=kerneldensity(b1glr[:,2])

### rule of thumb
function rtbwd(svec,n)
  return std(svec)*1.06*n^(-0.2)
end


l=@layout [a b; c d]

p1= plot(b1glr[:,4:6],[kerneldensity(b1glr[:,4],h=rtbwd(b1glr[:,4],n)) kerneldensity(b1glr[:,5],h=rtbwd(b1glr[:,5],n)) kerneldensity(b1glr[:,6],h=rtbwd(b1glr[:,6],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(a)")

p2= plot(b1lft3[:,4:6],[kerneldensity(b1lft3[:,4],h=rtbwd(b1lft3[:,4],n)) kerneldensity(b1lft3[:,5],h=rtbwd(b1lft3[:,5],n)) kerneldensity(b1lft3[:,6],h=rtbwd(b1lft3[:,6],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(b)")

p3= plot(b1glr[:,[2, 5, 8]],[kerneldensity(b1glr[:,2],h=rtbwd(b1glr[:,2],n)) kerneldensity(b1glr[:,5],h=rtbwd(b1glr[:,5],n)) kerneldensity(b1glr[:,8],h=rtbwd(b1glr[:,8],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(c)")

p4= plot(b1lft3[:,[2,5,8]],[kerneldensity(b1lft3[:,2],h=rtbwd(b1lft3[:,2],n)) kerneldensity(b1lft3[:,5],h=rtbwd(b1lft3[:,5],n)) kerneldensity(b1lft3[:,8],h=rtbwd(b1lft3[:,8],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(d)")

plot(p1,p2,p3,p4, layout=l)

savefig("p1a3.pdf")


### simulations for error distributions
### Figure 2

n= 200
nsim=1000
rind=2
err=[1 2 3 4]
#

b1glr=zeros(nsim,length(err))
#
b1lft1=zeros(nsim,length(err))
b1lft2=zeros(nsim,length(err))
b1lft3=zeros(nsim,length(err))
b1lft4=zeros(nsim,length(err))

# iterating over all the parameters
incr=1
#
for i in 1:length(err)

data=randomdata(nsim,n,0.6,1,err[i])
y=data[5][:,2]
xmat=hcat(data[1][:,2],data[2][:,2])
bw=[data[3][2],data[4][2]]
intgridmat=hcat(data[6][:,2],data[7][:,2])
#
# unrestr model
  uresm=optbwd(xmat,y,intgridmat)
  bwur=uresm["bwd"]

for j in 1:nsim

y=data[5][:,j]
xmat=hcat(data[1][:,j],data[2][:,j])
intgridmat=hcat(data[6][:,j],data[7][:,j])
#
# unrestr model
  uresm=ssb(xmat,y,intgridmat,bwur)
#  bwur=uresm["bwd"]
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


Gmd= resm["Gmat"]  #xmat_res*inv(xmat_res'*xmat_res)*xmat_res'
Gd= uresm["Gjmat"][:,:,rind]
Igmd=(I(n)-Gmd)*Gd
Pgmd=Igmd*inv(Igmd'*Igmd)*Igmd'
Ig=I(n)-uresm["Gmat"]
C= Gd+ Ig*uresm["Hstar"][:,:,rind]+uresm["Hstar"][:,:,rind]*Ig-Ig*uresm["Hstar"][:,:,rind]*uresm["Hstar"][:,:,rind]'*Ig #Pgmd
#D=Ig-reduce(+, [uresm["Hstar"][:,:,i]*Ig+Ig*uresm["Hstar"][:,:,i]- Ig*uresm["Hstar"][:,:,i]*uresm["Hstar"][:,:,i]*Ig for i in 1:nvar],dims=1)[1]
E=Gd+ Ig*uresm["Hstar"][:,:,rind]+uresm["Hstar"][:,:,rind]*Ig
#Fglr=((y'*C*y)/(y'*D*y))*(tr(D)/tr(C  Pgmd
EtE=E'*E

sn=0
dn=0
for i in 1:(n-1)
for j in (i+1):n
sn=sn+ C[i,j]^2
dn=dn+EtE[i,j]^2
end
end

# normalized lambda
rk=tr(C)/(2*sn)
# normalized q
sk=2*tr(EtE)/(2*dn)

#rkglr[j,incr]=rk
#sklft[j,incr]=sk


# storing normalized lambda and q statistics
b1glr[j,i]=lambda  #rk*
b1lft1[j,i]=Q1   #sk*
b1lft2[j,i]=Q2
b1lft3[j,i]=Q3
b1lft4[j,i]=Q4
# for j
end

# for i
end

##
using KernelEstimator
using Plots

td=kerneldensity(b1glr[:,2])

l=@layout [a b]
p1= plot(b1glr[:,1:4],[kerneldensity(b1glr[:,1],h=rtbwd(b1glr[:,1],n)) kerneldensity(b1glr[:,2],h=rtbwd(b1glr[:,2],n)) kerneldensity(b1glr[:,3],h=rtbwd(b1glr[:,3],n)) kerneldensity(b1glr[:,4],h=rtbwd(b1glr[:,4],n)) ],seriestype=:line,legend=false,linestyle=[ :solid :dash  :dot :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(a)")

p2= plot(0.5 .* b1lft1[:,1:4],[kerneldensity(0.5 .* b1lft1[:,1],h=rtbwd(0.5 .*b1lft1[:,1],n)) kerneldensity(b1lft1[:,2], h=rtbwd(b1lft1[:,2],n)) kerneldensity(b1lft1[:,3],h=rtbwd(b1lft1[:,3],n)) kerneldensity(b1lft1[:,4],h=rtbwd(b1lft1[:,4],n))  ],seriestype=:line,legend=false,linestyle=[:solid :dash  :dot :dashdot],xlims=(0, 100),ylims=(0, 0.07),title="(b)")

#
plot(p1,p2, layout=l)

savefig("p1b3.pdf")
