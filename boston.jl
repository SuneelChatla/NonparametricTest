#
using  Statistics, StatsBase, Plots, Random

#
include("ssb-inf.jl")

bdata=housing()



# regular matrix format: 1-13 are independent variables, 14 is dependent variable.
#
#    1. CRIM      per capita crime rate by town
#    2. ZN        proportion of residential land zoned for lots over
#                 25,000 sq.ft.
#    3. INDUS     proportion of non-retail business acres per town
#    4. CHAS      Charles River dummy variable (= 1 if tract bounds
#                 river; 0 otherwise)
#    5. NOX       nitric oxides concentration (parts per 10 million)
#    6. RM        average number of rooms per dwelling
#    7. AGE       proportion of owner-occupied units built prior to 1940
#    8. DIS       weighted distances to five Boston employment centres
#    9. RAD       index of accessibility to radial highways
#    10. TAX      full-value property-tax rate per $10,000
#    11. PTRATIO  pupil-teacher ratio by town
#    12. B        1000(Bk - 0.63)^2 where Bk is the proportion of blacks
#                 by town
#    13. LSTAT    % lower status of the population
#    14. MEDV     Median value of owner-occupied homes in $1000's
###

xfull=bdata[1]'
yfull=bdata[2]'

# selecting the same 4 variables from Fan & Jiang (2005)
xmat=hcat(xfull[:, 6],log.(xfull[:,10]),xfull[:,11], log.(xfull[:,13]))
y=yfull

 n, p=size(xmat)

Ngrid=401
intgridmat=zeros(Ngrid,p)
for j in 1:p
intgridmat[:,j]=range(minimum(xmat[:,j]),maximum(xmat[:,j]),length=Ngrid)
end

#
optmodel=optbwd(xmat,y,intgridmat)

# removing outliers
te=y .- optmodel["fit"]
noutind=findall(x -> (-11 < x < 12), vec(te))
#
xfull=xfull[noutind,:]
yfull=yfull[noutind]
### random sample of 200 observations
# Random.seed!(202103)
# ind=sample(1:size(xfull)[1],200)
# selecting the same 4 variables from Fan & Jiang (2005)
xmat=hcat(xfull[:, 6],log.(xfull[:,10]),xfull[:,11], log.(xfull[:,13]))
y=yfull

n, p=size(xmat)

 Ngrid=401
intgridmat=zeros(Ngrid,p)
for j in 1:p
intgridmat[:,j]=range(minimum(xmat[:,j]),maximum(xmat[:,j]),length=Ngrid)
end

#
optmodel=optbwd(xmat,y,intgridmat)

#
fjmodel=kcb(xmat,y,optmodel["bwd"])
#
#df=DataFrame(xmat,["x$i" for i in 1:size(xmat)[2]])
#  df[!,:y]=vec(y) .-mean(vec(y))
#h=optmodel["bwd"]
 # R"library(wsbackfit)"
 # R"m0=sback(formula = y ~ sb(x1, h = $h[1]) + sb(x2, h = $h[2]) + sb(x3, h=$h[3])+ sb(x4, h=$h[4]), $df)"
 # R"effects=m0$effects"
 # R"lincoef=m0$coeff"

 cxmat=(xmat .- mean(xmat,dims=1))
sbmodel=SnM2021(cxmat,y .- mean(y),optmodel["bwd"],2) # rind=2 ; optional



# using mgcv

rdf=DataFrame(y=vec(y),x1=xmat[:,1],x2=xmat[:,2],x3=xmat[:,3],x4=xmat[:,4])
using RCall
R"library(mgcv)"
R"gmod=gam(y~s(x1)+s(x2)+s(x3)+s(x4),data=$rdf)"
R"str(gmod)"
R"""
gmodplot=predict(gmod,type="terms")
"""
@rget gmodplot

R"cf=coef(gmod)"
@rget cf
##

# testplot
using Plots,LaTeXStrings
plot(xmat,optmodel["mu"][:,2:5],seriestype=:line, layout=(2,2))
plot(xmat,optmodel["parmet"][:,2:5],seriestype=:line, layout=(2,2))

##
# partial residual plots
parred=zeros(n,p)
parred_fj=zeros(n,p)
parred_mg=zeros(n,p)
#
for j in 1:p
parred[:,j]= y .- sum(optmodel["mu"][:,1:end .!=j+1],dims=2)
parred_fj[:,j]=  y .-fill(mean(y),n) - sum(fjmodel["mu"][:,1:end .!=j+1],dims=2)
parred_mg[:,j]= y .- fill(cf[1],n) - sum(gmodplot[:,1:end .!=j+1],dims=2)
end

# plot for ssbf
mu= optmodel["mu"]
pbeta=optmodel["parmet"]
mustar=optmodel["mustar"]

l=@layout [a b; c d]

plot(xmat[:,1],parred[:,1],seriestype=:scatter,legend=:bottomright,xlabel="RM",ylabel="additive residual",title="",ylims=[-15,20],seriesalpha=0.15,label="");

plot!(xmat[:,1], mu[:,2],seriestype=:line,legend=:bottomright,linestyle=:solid,lw=2,label=L"\widehat{m}");

plot!(xmat[:,1], pbeta[:,2],seriestype=:line,legend=:bottomright,linestyle=:dash,lw=2,label=L"\widehat{g}");

p1=plot!(xmat[:,1], mustar[:,1],seriestype=:line,legend=:bottomright,linestyle=:dot,lw=2,label=L"\widehat{m}^*");
#

plot(xmat[:,2],parred[:,2],seriestype= :scatter,legend=:bottomright,xlabel="log(TAX)",ylabel="additive residual",label="",seriesalpha=0.15);

plot!(xmat[:,2], mu[:,3],seriestype=:line,legend=:bottomright,linestyle=:solid,lw=2,label=L"\widehat{m}");

plot!(xmat[:,2], pbeta[:,3],seriestype=:line,legend=:bottomright,linestyle=:dash,lw=2,label=L"\widehat{g}");

p2=plot!(xmat[:,2], mustar[:,2],seriestype=:line,legend=:bottomright,linestyle=:dot,lw=2,label=L"\widehat{m}^*");
#

plot(xmat[:,3],parred[:,3],seriestype=:scatter,legend=:bottomright,xlabel="PTRATIO",ylabel="additive residual",label="",seriesalpha=0.15);

plot!(xmat[:,3], mu[:,4],seriestype=:line,legend=:bottomright,linestyle=:solid,lw=2,label=L"\widehat{m}");

plot!(xmat[:,3], pbeta[:,4],seriestype=:line,legend=:bottomright,linestyle=:dash,lw=2,label=L"\widehat{g}");

p3=plot!(xmat[:,3], mustar[:,3],seriestype=:line,legend=:bottomright,linestyle=:dot,lw=2,label=L"\widehat{m}^*");
#
plot(xmat[:,4], parred[:,4] ,seriestype=:scatter,legend=:topright,xlabel="log(LSTAT)",ylabel="additive residual",label="",seriesalpha=0.15);

plot!(xmat[:,4], mu[:,5],seriestype=:line,legend=:topright,linestyle=:solid,lw=2,label=L"\widehat{g}");

plot!(xmat[:,4], pbeta[:,5],seriestype=:line,legend=:topright,linestyle=:dash,lw=2,label=L"\widehat{g}");

p4=plot!(xmat[:,4], mustar[:,4],seriestype=:line,legend=:topright,linestyle=:dot,lw=2,label=L"\widehat{m}^*");


plot(p1,p2,p3,p4, layout=l)
plot!(size=(950,750))
savefig("b500pr.pdf")


#
# plot for ssbf
mu_fj= fjmodel["mu"]
mu_mg=gmodplot
bfsn, npcoef, pcoef=sbmodel
mu_sb=  cxmat .* pcoef[2:length(pcoef)]' .+ npcoef

l=@layout [a b; c d]

plot(xmat[:,1],parred[:,1],seriestype=:scatter,legend=:bottomright,xlabel="RM",ylabel="additive residual",seriesalpha=0.15,ylims=[-15,20],seriescolor=:gray,label="");

plot!(xmat[:,1], mu[:,2],seriestype=:line,legend=:bottomright,linestyle=:solid,lw=2,label="SSB");

plot!(xmat[:,1], mu_fj[:,1],seriestype=:line,legend=:bottomright,linestyle=:dash,lw=2,label="CB");

p11=plot!(xmat[:,1], mu_mg[:,1],seriestype=:line,legend=:bottomright,linestyle=:dot,lw=2 ,label="Penalized");

plot!(xmat[:,1], mu_sb[:,1],seriestype=:line,legend=:bottomright,linestyle=:dashdot,lw=2,label="SB")
#

plot(xmat[:,2],parred[:,2],seriestype= :scatter,legend=:bottomright,xlabel="log(TAX)",ylabel="additive residual",seriesalpha=0.15,ylims=[-15,20],seriescolor=:gray,label="");

plot!(xmat[:,2], mu[:,3],seriestype=:line,legend=:bottomright,linestyle=:solid,lw=2,label="SSB");

plot!(xmat[:,2], mu_fj[:,2],seriestype=:line,legend=:bottomright,linestyle=:dash,lw=2,label="CB");

p21=plot!(xmat[:,2], mu_mg[:,2],seriestype=:line,legend=:bottomright,linestyle=:dot,lw=2,label="Penalized");

plot!(xmat[:,2], mu_sb[:,2],seriestype=:line,legend=:bottomright,linestyle=:dashdot,lw=2,label="SB")
#,label="SB"
#

plot(xmat[:,3],parred[:,3],seriestype=:scatter,legend=:bottomright,xlabel="PTRATIO",ylabel="additive residual",title="",seriesalpha=0.15,ylims=[-15,20],seriescolor=:gray,label="");

plot!(xmat[:,3], mu[:,4],seriestype=:line,legend=:bottomright,linestyle=:solid,lw=2,label="SSB");

plot!(xmat[:,3], mu_fj[:,3],seriestype=:line,legend=:bottomright,linestyle=:dash,lw=2,label="CB");

p31=plot!(xmat[:,3], mu_mg[:,3],seriestype=:line,legend=:bottomright,linestyle=:dot,lw=2,label="Penalized");

plot!(xmat[:,3], mu_sb[:,3],seriestype=:line,legend=:bottomright,linestyle=:dashdot,lw=2,label="SB")

#
plot(xmat[:,4], parred[:,4] ,seriestype=:scatter,legend=:topright,xlabel="log(LSTAT)",ylabel="additive residual",seriesalpha=0.15,ylims=[-15,20],seriescolor=:gray,label="");

plot!(xmat[:,4], mu[:,5],seriestype=:line,legend=:topright,linestyle=:solid,lw=2,label="SSB");

plot!(xmat[:,4], mu_fj[:,4],seriestype=:line,legend=:topright,linestyle=:dash,lw=2,label="CB");

p41=plot!(xmat[:,4], mu_mg[:,4],seriestype=:line,legend=:topright,linestyle=:dot,lw=2,label="Penalized");

plot!(xmat[:,4], mu_sb[:,4],seriestype=:line,legend=:topright,linestyle=:dashdot,lw=2,label="SB")

plot(p11,p21,p31,p41, layout=l)
plot!(size=(900,600))
savefig("b500pr_cbf_mgc_sb.pdf")


## conditional bootstrap for power analysis
### for hypothesis 1 (same as in Fan 2005)
# different bandwidth choices
ncboot=1000
bwopt=optmodel["bwd"]
bwdcom=zeros(5,p)

bwdcom[1,:]=bwopt ./ 2
bwdcom[2,:]=2 .* bwopt ./3
bwdcom[3,:]=bwopt
bwdcom[4,:]=3 .* bwopt ./ 2
bwdcom[5,:]=2 .* bwopt
#
rss1=zeros(5)
rss0=zeros(5)
rss0_fj=zeros(5)
rss1_fj=zeros(5)
#
lambda=zeros(5)
Q1=zeros(5)
Q2=zeros(5)
Q3=zeros(5)
Q4=zeros(5)
Flambda=zeros(5)
Fq=zeros(5)
lambda_fj=zeros(5)
sn=zeros(5)
#
bplambda=zeros(5)
bpq1=zeros(5)
bpq2=zeros(5)
bpq3=zeros(5)
bpq4=zeros(5)
bplambda_fj=zeros(5)
bpsn=zeros(5)
bpFlambda=zeros(5)
bpFq=zeros(5)
#
lambdastar=zeros(ncboot,5)
Q1star= zeros(ncboot,5)
Q2star= zeros(ncboot,5)
Q3star= zeros(ncboot,5)
Q4star= zeros(ncboot,5)
lambdastar_fj=zeros(ncboot,5)
snstar=zeros(ncboot,5)
#
Flambdastar=zeros(ncboot,5)
Fqstar=zeros(ncboot,5)

# loop for bandwidths
j=1
#for j in 1:5
        # calculation of test statistics for the original data
        bwur=bwdcom[j,:]
        urmodel=ssb(xmat,y,intgridmat,bwur)
        fitur=urmodel["fit"]
        #bwur=uresm["bwd"]
        # restricted model data
        rind=[2 3 4]
        xmat_res=xmat[:, 1]
        intgridmat_res=intgridmat[:,1]
        # h star matrix for the nonparametric component 1
        Hmat1=hstar(xmat_res,intgridmat_res,bwur[1])[1]
        # for parametric null hypothesis
        xcmat_rp=reshape(xmat[:,rind],n,length(rind))
        xmat_res_pm=xcmat_rp .- mean(xcmat_rp,dims=1) #hcat(ones(n,1), xcmat_rp .- mean(xcmat_rp,dims=1))
        # two step semiparametric model
        betahat=xmat_res_pm \ (y .- mean(y))
        mu1=Hmat1*(y-xmat_res_pm*betahat .- mean(y))
        betahat=xmat_res_pm \(y-mu1 .- mean(y))
        mu1=Hmat1*(y-xmat_res_pm*betahat .- mean(y))
        #

        fitr=mu1+xmat_res_pm*betahat .+ mean(y)

        rss1[j]=sum((y .- fitur) .^ 2)
        rss0[j]=sum((y .- fitr) .^ 2)

        ### chi-square statistics
        lambda[j]= n/2 * (rss0[j]-rss1[j])/rss1[j]
        Q1[j]= n* sum(linex.(1e-6,1,fitur-fitr))/rss1[j]
        Q2[j]= n* sum(linex.(0.2,1,fitur-fitr))/rss1[j]
        Q3[j]= n* sum(linex.(0.5,1,fitur-fitr))/rss1[j]
        Q4[j]= n* sum(linex.(1,1,fitur-fitr))/rss1[j]

        # F test statistic calculation
        xGm1=hcat(ones(n),xcmat_rp)
        Gm1=xGm1*inv(xGm1'*xGm1)*xGm1'
        IGm1=I(n)-Gm1
        PIgm11=IGm1*xmat_res*inv(xmat_res'*IGm1*IGm1*xmat_res)*xmat_res'*IGm1
        xxmat=hcat(ones(n),xmat)
        Gmat=xxmat*inv(xxmat'*xxmat)*xxmat'
        Hresmat=sum(urmodel["Hstar"][:,:,2:4],dims=3)

        C=PIgm11 .+ (I(n) .- Gmat)* Hresmat[:,:,1] .+ Hresmat[:,:,1]*(I(n) .- Gmat)-reduce(+, [(I(n) .- Gmat)* urmodel["Hstar"][:,:,i] *urmodel["Hstar"][:,:,i]*(I(n) .- Gmat) for i in 2:size(xmat)[2]])
        E=PIgm11 .+ Hresmat[:,:,1]*(I(n) .- Gmat)
        D=(I(n) .- Gmat)-reduce(+, [urmodel["Hstar"][:,:,i]*(I(n) .- Gmat) for i in 1:size(xmat)[2]],dims=1)[1]
        #
        Flambda[j]= (y'*C*y)/(y'*D*y)*(tr(D)/tr(C))
        Fq[j]=(y'*E'*E*y)/(y'*D*y)*(tr(D)/tr(E'*E))
        #
        epsi_ures=y .- fitur
        cen_epsi_ures=epsi_ures .- mean(epsi_ures)

        ## fj model setup
        urmodel_fj=kcb(xmat,y,bwur)
        fitur_fj=urmodel_fj["fit"]
        #bwur=uresm["bwd"]
        # S matrix for the nonparametric component 1
        Smat1=Smat(xmat_res,xmat_res,bwur[1])[1]
        # for parametric null hypothesis
        #xcmat_rp=reshape(xmat[:,rind],n,length(rind))
        #xmat_res_pm=hcat(ones(n,1), xcmat_rp .- mean(xcmat_rp,dims=1))
        # two step semiparametric model
        betahat=xmat_res_pm \ y
        mu1_fj=Smat1*(y-xmat_res_pm*betahat .- mean(y))
        betahat=xmat_res_pm \(y-mu1_fj .- mean(y))
        mu1_fj=Smat1*(y-xmat_res_pm*betahat .- mean(y))
        #

        fitr_fj=mu1_fj+xmat_res_pm*betahat .+ mean(y)

        rss1_fj[j]=sum((y .- fitur_fj) .^ 2)
        rss0_fj[j]=sum((y .- fitr_fj) .^ 2)

        ### chi-square statistics
        lambda_fj[j]= n/2 * (rss0_fj[j]-rss1_fj[j])/rss1_fj[j]

        epsi_ures_fj=y .- fitur_fj
        cen_epsi_ures_fj=epsi_ures_fj .- mean(epsi_ures_fj)

        ########### 
        ### Sn statistic 
        h=bwur
            # creatiin of dataframe
        xcmat=xmat .- mean(xmat,dims=1)
        df=DataFrame(xcmat,["x$i" for i in 1:size(xmat)[2]])
        df[!,:y]=vec(y) .- mean(y)
        #  using wsbackfit package to do the calcualtions
        R"library(wsbackfit)"
        R"m0=sback(formula = y ~ sb(x1, h = $h[1]) + sb(x2, h = $h[2]) + sb(x3, h=$h[3])+ sb(x4, h=$h[4]), $df)"
        R"effects=m0$effects"
        R"lincoef=m0$coeff"
        
        # return to julia
        @rget(effects)
        @rget(lincoef)
        # computation of combined effect (parametric + nonparametric)
        m=effects #+ xcmat .* lincoef[2:length(lincoef)]'

        snvec=zeros(length(rind))
        # di=1
        for di in 1:length(rind)
            # density weights calcuation
            p2hat=kerneldensity(xmat[:,rind[di]],h=h[rind[di]])
            # sorting
            tmat=hcat(xmat[:,rind[di]],m[:,rind[di]],p2hat)
            tomat=tmat[sortperm(tmat[:,1]),:]
            # Riemann weights
            diffvec=zeros(length(xmat[:,rind[di]]))
            diffvec[2:size(xmat)[1]]=[abs(tomat[i-1,1]-tomat[i,1])  for i in 2:size(xmat)[1]]
            # computing the statistic
            snvec[di]=sum((tomat[:,2] .^2) .* tomat[:,3] .* diffvec) 
        end
        #
        sn[j] = sum(snvec)

        fitr_sn=rcopy(R"sback(formula = y ~ sb(x1, h = $h[1]) + x2+ x3+ x4, $df)$fitted.values") .+ mean(y)

        epsi_ures_sn=rcopy(R"m0$residuals")
        cen_epsi_ures_sn=epsi_ures_sn .- mean(epsi_ures_sn)

        #
        ################################
        # loop through number of bootstraps
        #
        #i=1
        for i in 1:100 # 1:ncboot
            ystar =  fitr .+ sample(cen_epsi_ures,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
            # bootstrapped
            uresmstar= ssb(xmat,ystar,intgridmat,bwur)
            #
            betahatstar=xmat_res_pm \ (ystar .- mean(ystar))
            mu1star=Hmat1*(ystar-xmat_res_pm*betahatstar .- mean(ystar))
            betahatstar=xmat_res_pm \(ystar-mu1star .- mean(ystar))
            mu1star=Hmat1*(ystar-xmat_res_pm*betahatstar .- mean(ystar))

            #
            fitrstar=mu1star+xmat_res_pm*betahatstar .+ mean(ystar)
            rss1star=sum((ystar .- uresmstar["fit"]) .^ 2)
            rss0star=sum((ystar .- fitrstar) .^ 2)

            lambdastar[i,j]= n/2 * (rss0star-rss1star)/rss1star
            #
            Q1star[i,j]= n* sum(linex.(1e-6,1,uresmstar["fit"]-fitrstar))/rss1star
            Q2star[i,j]= n* sum(linex.(0.2,1,uresmstar["fit"]-fitrstar))/rss1star
            Q3star[i,j]= n* sum(linex.(0.5,1,uresmstar["fit"]-fitrstar))/rss1star
            Q4star[i,j]= n* sum(linex.(1,1,uresmstar["fit"]-fitrstar))/rss1star
            ### Ftest calculation
            
            Flambdastar[i,j]= (ystar'*C*ystar)/(ystar'*D*ystar)*(tr(D)/tr(C))
            Fqstar[i,j]=(ystar'*E'*E*ystar)/(ystar'*D*ystar)*(tr(D)/tr(E'*E))
            
            ###
            ystar_fj =  fitr_fj .+ sample(cen_epsi_ures_fj,n) #sum(uresm["mu"][:,1:end .!= rind+1],dims=2)
            # bootstrapped
            uresmstar_fj= kcb(xmat,ystar_fj,bwur)
            
            betahatstar=xmat_res_pm \ (ystar_fj .- mean(ystar_fj))
            mu1star_fj=Smat1*(ystar_fj-xmat_res_pm*betahatstar .- mean(ystar_fj))
            betahatstar=xmat_res_pm \(ystar_fj-mu1star_fj .- mean(ystar_fj))
            mu1star_fj=Smat1*(ystar_fj-xmat_res_pm*betahatstar .- mean(ystar_fj))

            #
            fitrstar_fj=mu1star_fj+xmat_res_pm*betahatstar .+ mean(ystar_fj)
            rss1star_fj=sum((ystar_fj .- uresmstar_fj["fit"]) .^ 2)
            rss0star_fj=sum((ystar_fj .- fitrstar_fj) .^ 2)

            lambdastar_fj[i,j]= n/2 * (rss0star_fj-rss1star_fj)/rss1star_fj

            ###### Sn test statistics
            ystar_sn =  fitr_sn .+ sample(cen_epsi_ures_sn,n)

            df[!,:ystar]=vec(ystar_sn) .- mean(ystar_sn)
            #  using wsbackfit package to do the calcualtions
            R"library(wsbackfit)"
            R"m0=sback(formula = ystar ~ sb(x1, h = $h[1]) + sb(x2, h = $h[2]) + sb(x3, h=$h[3])+ sb(x4, h=$h[4]), $df)"
            R"effects=m0$effects"
            R"lincoef=m0$coeff"
        
        # return to julia
            @rget(effects)
            @rget(lincoef)
            # computation of combined effect (parametric + nonparametric)
            m=effects # + xcmat .* lincoef[2:length(lincoef)]'

            snvec=zeros(length(rind))
            # di=1
            for di in 1:length(rind)
                # density weights calcuation
                p2hat=kerneldensity(xmat[:,rind[di]],h=h[rind[di]])
                # sorting
                tmat=hcat(xmat[:,rind[di]],m[:,rind[di]],p2hat)
                tomat=tmat[sortperm(tmat[:,1]),:]
                # Riemann weights
                diffvec=zeros(length(xmat[:,rind[di]]))
                diffvec[2:size(xmat)[1]]=[abs(tomat[i-1,1]-tomat[i,1])  for i in 2:size(xmat)[1]]
                # computing the statistic
                snvec[di]=sum((tomat[:,2] .^2) .* tomat[:,3] .* diffvec) 
            end
            #
            snstar[i,j] = sum(snvec)
        end
        #
        using Tables
        cbstat=Tables.table(hcat(lambdastar,Q1star, Q2star, Q3star, Q4star,Flambdastar,Fqstar,lambdastar_fj,snstar))
        #
        bplambda[j]=mean(lambda[j] .<= lambdastar[:,j])
        bpq1[j]=mean(Q1[j] .<= Q1star[:,j])
        bpq2[j]=mean(Q2[j] .<= Q2star[:,j])
        bpq3[j]=mean(Q3[j] .<= Q3star[:,j])
        bpq4[j]=mean(Q4[j] .<= Q4star[:,j])

        bpFlambda[j] = mean(Flambda[j] .<= Flambdastar)
        bpFq[j] = mean(Fq[j] .<= Fqstar)

        bplambda_fj[j]=mean(lambda_fj[j] .<= lambdastar_fj[:,j])
        #
        bpsn[j] = mean(sn[j] .<= snstar[:,j])
        # end of bootstrap
end

##