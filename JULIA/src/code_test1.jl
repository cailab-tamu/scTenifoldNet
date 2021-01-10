cd(dirname(@__FILE__))

include("ScTenifoldNet.jl")
using .ScTenifoldNet
X0=rand(100,1000);
X1=copy(X0)
X1[4,:].=0.0


d,fc,p,adjp=ScTenifoldNet.sctenifoldnet(X0,X1,donorm=false)

#@show Threads.nthreads()
#@time Z0=ScTenifoldNet.tenrnet(X0, donorm=false)
#@time Z1=ScTenifoldNet.tenrnet(X1, donorm=false)
#Z0=0.5*(Z0+Z0')
#Z1=0.5*(Z1+Z1')
#@time d,aln0,aln1=ScTenifoldNet.manialn(Z0,Z1)
#fc,p,adjp=ScTenifoldNet.drgenes(d)

#using StatsPlots, Distributions
#x=rand(Chisq(1), length(fc)) 
#qqplot(x, fc)