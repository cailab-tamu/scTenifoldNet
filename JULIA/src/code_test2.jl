cd(dirname(@__FILE__))

include("ScTenifoldNet.jl")
using .ScTenifoldNet

using LinearAlgebra, Statistics, Distributions, Random

d=NegativeBinomial(20,0.98)
X=rand(d,100,2000)
lbszv=30

Y=copy(X)
Y[10,:]=Y[50,:]
Y[2,:]=Y[11,:]
Y[3,:]=Y[5,:]

X=X[:,vec(sum(X,dims=1).>lbszv)]
Y=Y[:,vec(sum(Y,dims=1).>lbszv)]

@show Threads.nthreads()

@time Z0=ScTenifoldNet.tenrnet(X, donorm=true)
@time Z1=ScTenifoldNet.tenrnet(Y, donorm=true)
Z0=0.5*(Z0+Z0')
Z1=0.5*(Z1+Z1')
@time d=ScTenifoldNet.manialn(Z0,Z1)
fc,p,adjp=ScTenifoldNet.drgenes(d)

@show [adjp[10] adjp[11] adjp[50] adjp[20]];
@show [adjp[70] adjp[31] adjp[55] adjp[26]];

using StatsPlots
x=rand(Chisq(1), length(fc)) 
qqplot(x, fc)
