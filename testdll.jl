using Statistics, StatsBase, LinearAlgebra
############## test for varying coeff estimation: sbvcest ##########
# set Integers variables
n=150       # Number of observations
d=3         # Number of dimensions
ng=40       # Grid length (not too small, used for integration)
maxit=30    # Maximal number of iterations
m=2         # for trimming in test, see above
nb=500      # number of bootstrap samples
# set Matrix/Array/Vector
x=rand(n,d)     # Design matrix X for nonparametric components
mf= zeros(ng,d)               # estimates on grid under H1 (output)
z=randn(n,d) .+ 2.0   # Design matrix Z, linear multipliers
h=Diagonal([std(x[:,i]) for i=1:3] ./ n^(1/5) .* 1.5) # Bandwidth vector for estimation
# set Double precision
conv=0.0001   # Convergence criteria in backfitting
hd=std(x[:,1])/n^(1/5)   # bandwidth for test statistic
to=999        # value of test statistic
m0=0          # output-constant of model (c in notation below)
rd=1000.0    # random seed for replicability of simulations
pv=1.0        # initialize p-value 
# simulate a varying coefficient model Y = sum_j b_j(x_j) z_j + error
y=zeros(n)
 y = (sin.(pi .* x[:,2]) .+ 2) .* z[:,2] .+ (sin.(pi .* x[:,3]) .+ 2) .* z[:,3] .+ randn(n) 

 xg = mf; 
 for j in 1:d  # equispaced grid 
   s = ( maximum(x[:,j]) .- minimum(x[:,j])) ./ (ng-1)
   for i in 1:ng
     xg[i,j] = minimum(x[:,j]) .+ (i - 1) .* s 
   end
end   
# call first routine, assign to list erg
erg<-.Fortran("sbvctest",x=x,z=z,y=y,n=n,d=d,h=h,mf=mf,xg=xg
           ,ng=ng,maxit=maxit,conv=conv,m=m,hd=hd,to=to,m0=m0)