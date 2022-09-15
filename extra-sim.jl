
include("ssb-inf.jl")

##
#### test code
### Assume equally-spaced grid points  of size 401 on [-1,1] for either domemsion for estimation
nobs=[200 500]

ed=[1 2 3 4]
theta=[0 0.2 0.4 0.6 0.8 1.0]
ncboot=3000
nsim=1000
#
incr=1
onemat=zeros(length(nobs)*length(ed)*length(theta),6)
fivemat=zeros(length(nobs)*length(ed)*length(theta),6)
tenmat=zeros(length(nobs)*length(ed)*length(theta),6)

# iterating over all the parameters
for i in 1:length(nobs)
bwur=0.45 .* fill(nobs[i]^(-1/5),4)
for j in 1:length(ed)
for k in 1:length(theta)

data=randomdatarev(nsim,nobs[i],0.6,theta[k],ed[i])
#

#bwur= [0.15 0.15] #optbwd(xmat,y,intgridmat)["bwd"]
res=cbootrev(data,ncboot,bwur,2)
#
onemat[incr,:]=res[1]
fivemat[incr,:]=res[2]
tenmat[incr,:]=res[3]
global incr=incr+1
#
#end
#end
end
# write the results to csv
CSV.write("ones5.csv",convert(DataFrame,onemat))
CSV.write("fives5.csv",convert(DataFrame,fivemat))
CSV.write("tens5.csv",convert(DataFrame,tenmat))

####
##  Figure 1
### Wilks phenomenon

include("ssb-inf.jl")

#
n= 200

beta=[-1.5 0 1.5]
nsim=1000
rind=2
bwdcom=zeros(3,4)
#

b1glr=zeros(nsim,length(beta)*3)
rkglr=zeros(nsim,length(beta)*3)
b3glr=zeros(nsim,length(beta)*3)
#
b1lft1=zeros(nsim,length(beta)*3)
b1lft2=zeros(nsim,length(beta)*3)
b1lft3=zeros(nsim,length(beta)*3)
b1lft4=zeros(nsim,length(beta)*3)
#
b1glr_fj=zeros(nsim,length(beta)*3)


sklft=zeros(nsim,length(beta)*3)


# iterating over all the parameters
incr=1
#
for i in 1:length(beta)

data=randomdatarev(nsim,n,0.6,1,1,beta[i])
y=data[9][:,2]
xmat=hcat(data[1][:,2],data[2][:,2],data[3][:,2],data[4][:,2])
bw=[data[5][2],data[6][2],data[7][2],data[8][2]]
intgridmat=hcat(data[10][:,2],data[11][:,2],data[12][:,2],data[13][:,2])
#
# unrestr model
uresm=optbwd(xmat,y,intgridmat)
bw=uresm["bwd"]
bwdcom[2,:]=bw'
bwdcom[1,:]=[bw[1]*1/3,bw[2],bw[3],bw[4]]
bwdcom[3,:]=[bw[1]*1.5,bw[2],bw[3],bw[4]]


#

for k in 1:3
#
bwur=bwdcom[k,:]
Threads.@threads for j in 1:nsim

y=data[9][:,j]
xmat=hcat(data[1][:,j],data[2][:,j],data[3][:,j],data[4][:,j])
intgridmat=hcat(data[10][:,j],data[11][:,j],data[12][:,j],data[13][:,j])
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

   # for fj
   uresm_fj=kcb(xmat,y,bwur)
#  bwur=uresm["bwd"]
  fitur_fj=uresm_fj["fit"]
  #bwur=uresm["bwd"]
  # restricted model data
  resm_fj=kcb(xmat_res,y,bw_res)
  fitr_fj=resm_fj["fit"]

  rss1_fj=sum((y .- fitur_fj) .^ 2)
  rss0_fj=sum((y .- fitr_fj) .^ 2)

# chi-square statistics
  lambda_fj= n/2 * (rss0_fj-rss1_fj)/rss1_fj



   #

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
sk=2*tr(EtE)/(2*dn)

rkglr[j,incr]=rk
sklft[j,incr]=sk


# storing normalized lambda and q statistics
b1glr[j,incr]=lambda  #rk*
b1lft1[j,incr]=Q1   #sk*
b1lft2[j,incr]=Q2
b1lft3[j,incr]=Q3
b1lft4[j,incr]=Q4
#
b1glr_fj[j,incr]=lambda_fj


# for j
end
# for k
incr=incr+1
end
# for i
end

##
CSV.write("wilks.csv", convert(DataFrame,hcat(b1glr,b1lft1,b1lft2,b1lft3,b1lft4,b1glr_fj)))
##
wilks=CSV.read("wilks.csv")

#td=kerneldensity(b1glr[:,2])
#
b1glr=Matrix(wilks[:,1:9])
b1lft3=Matrix(wilks[:,10:18])
b1glr_fj=Matrix(wilks[:,46:54])
##
using KernelEstimator
using   Plots

td=kerneldensity(b1glr[:,2])

### rule of thumb
using Statistics
function rtbwd(svec,n)
  return std(svec)*1.06*n^(-0.2)
end


l=@layout [a b c;  d e f]

p1= plot(b1glr[:,4:6],[kerneldensity(b1glr[:,4],h=rtbwd(b1glr[:,4],n)) kerneldensity(b1glr[:,5],h=rtbwd(b1glr[:,5],n)) kerneldensity(b1glr[:,6],h=rtbwd(b1glr[:,6],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(a)")

p2= plot(b1lft3[:,4:6],[kerneldensity(b1lft3[:,4],h=rtbwd(b1lft3[:,4],n)) kerneldensity(b1lft3[:,5],h=rtbwd(b1lft3[:,5],n)) kerneldensity(b1lft3[:,6],h=rtbwd(b1lft3[:,6],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(b)")



p3= plot(b1glr[:,[2, 5, 8]],[kerneldensity(b1glr[:,2],h=rtbwd(b1glr[:,2],n)) kerneldensity(b1glr[:,5],h=rtbwd(b1glr[:,5],n)) kerneldensity(b1glr[:,8],h=rtbwd(b1glr[:,8],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(d)")

p4= plot(b1lft3[:,[2,5,8]],[kerneldensity(b1lft3[:,2],h=rtbwd(b1lft3[:,2],n)) kerneldensity(b1lft3[:,5],h=rtbwd(b1lft3[:,5],n)) kerneldensity(b1lft3[:,8],h=rtbwd(b1lft3[:,8],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(e)")
#
p5= plot(b1glr_fj[:,4:6],[kerneldensity(b1glr_fj[:,4],h=rtbwd(b1glr_fj[:,4],n)) kerneldensity(b1glr_fj[:,5],h=rtbwd(b1glr_fj[:,5],n)) kerneldensity(b1glr_fj[:,6],h=rtbwd(b1glr_fj[:,6],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(c)")

p6= plot(b1glr_fj[:,[2, 5, 8]],[kerneldensity(b1glr_fj[:,2],h=rtbwd(b1glr_fj[:,2],n)) kerneldensity(b1glr_fj[:,5],h=rtbwd(b1glr_fj[:,5],n)) kerneldensity(b1glr_fj[:,8],h=rtbwd(b1glr_fj[:,8],n)) ],seriestype=:line,legend=false,linestyle=[:dash :solid :dashdot],xlims=(0, 60),ylims=(0, 0.07),title="(f)")




plot(p1,p2,p5,p3,p4,p6, layout=l)

using Plots
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
b1glr_fj=zeros(nsim,length(err))

# iterating over all the parameters
incr=1
#
Threads.@threads for i in 1:length(err)

data=randomdatarev(nsim,n,0.6,1,err[i])
y=data[9][:,2]
xmat=hcat(data[1][:,2],data[2][:,2],data[3][:,2],data[4][:,2])
bw=[data[5][2],data[6][2],data[7][2],data[8][2]]
intgridmat=hcat(data[10][:,2],data[11][:,2],data[12][:,2],data[13][:,2])
#
# unrestr model
  uresm=optbwd(xmat,y,intgridmat)
  bwur=uresm["bwd"]
#
for j in 1:nsim

y=data[9][:,j]
xmat=hcat(data[1][:,j],data[2][:,j],data[3][:,j],data[4][:,j])
intgridmat=hcat(data[10][:,j],data[11][:,j],data[12][:,j],data[13][:,j])
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

   # fj
   uresm_fj=kcb(xmat,y,bwur)
   #  bwur=uresm["bwd"]
     fitur_fj=uresm_fj["fit"]
     #bwur=uresm["bwd"]
     # restricted model data
     resm_fj=kcb(xmat_res,y,bw_res)
     fitr_fj=resm_fj["fit"]
   
     rss1_fj=sum((y .- fitur_fj) .^ 2)
     rss0_fj=sum((y .- fitr_fj) .^ 2)
   
   # chi-square statistics
     lambda_fj= n/2 * (rss0_fj-rss1_fj)/rss1_fj
   


   #

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
b1glr_fj[j,i]=lambda_fj
# for j
end

# for i
end

##
CSV.write("error.csv", convert(DataFrame,hcat(b1glr,b1lft1,b1lft2,b1lft3,b1lft4,b1glr_fj)))

##
using CSV,DataFrames

error=CSV.read("error.csv")

using KernelEstimator
using Plots

td=kerneldensity(b1glr[:,2])
#
b1glr=Matrix(error[:,1:4])
b1lft1=Matrix(error[:,5:8])
b1glr_fj=Matrix(error[:,21:24])
#
l=@layout [a b c]
p1= plot(b1glr[:,1:4],[kerneldensity(b1glr[:,1],h=rtbwd(b1glr[:,1],n)) kerneldensity(b1glr[:,2],h=rtbwd(b1glr[:,2],n)) kerneldensity(b1glr[:,3],h=rtbwd(b1glr[:,3],n)) kerneldensity(b1glr[:,4],h=rtbwd(b1glr[:,4],n)) ],seriestype=:line,legend=false,linestyle=[ :solid :dash  :dot :dashdot],xlims=(0, 70),ylims=(0, 0.07),title="(a)")


p2= plot( b1lft1[:,1:4],[kerneldensity( b1lft1[:,1],h=rtbwd(b1lft1[:,1],n)) kerneldensity(b1lft1[:,2], h=rtbwd(b1lft1[:,2],n)) kerneldensity(b1lft1[:,3],h=rtbwd(b1lft1[:,3],n)) kerneldensity(b1lft1[:,4],h=rtbwd(b1lft1[:,4],n))  ],seriestype=:line,legend=false,linestyle=[:solid :dash  :dot :dashdot],xlims=(0, 70),ylims=(0, 0.07),title="(b)")
#

p3= plot(b1glr_fj[:,1:4],[kerneldensity(b1glr_fj[:,1],h=rtbwd(b1glr_fj[:,1],n)) kerneldensity(b1glr_fj[:,2],h=rtbwd(b1glr_fj[:,2],n)) kerneldensity(b1glr_fj[:,3],h=rtbwd(b1glr[:,3],n)) kerneldensity(b1glr_fj[:,4],h=rtbwd(b1glr_fj[:,4],n)) ],seriestype=:line,legend=false,linestyle=[ :solid :dash  :dot :dashdot],xlims=(0, 70),ylims=(0, 0.07),title="(c)")
#
plot(p1,p2,p3, layout=l)

savefig("p1b3.pdf")




##########################################################################
######## Eigenvalues
##########################################################################

# data 
