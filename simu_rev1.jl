################################################################################
########  Data files are copied from Server
################################################################################
using CSV, DataFrames, CairoMakie
five= Matrix(CSV.read("oft5r1.csv",DataFrame))
sevent=Matrix(CSV.read("oft17rev1.csv",DataFrame))

####
#latex export
using LatexPrint

nobs=[100]

ed=[1 2 3 4]
theta=[0 0.2 0.4 0.6 0.8 1.0]

## oldcode
thetavec=repeat(theta',length(nobs)*length(ed))
edvec=repeat(vcat(repeat(ed,length(theta))...),1)
samplevec=vcat(repeat(nobs,24)...)

#
myfivemat=zeros(24,11)
myfivemat=hcat(edvec,thetavec,sevent[:,10:18])
tabular(myfivemat)

#
mytenmat=zeros(48,7)
mytenmat=hcat(edvec,thetavec,ten5)
tabular(mytenmat)

