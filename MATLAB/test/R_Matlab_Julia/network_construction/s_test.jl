using DelimitedFiles
cd(dirname(@__FILE__))

X=readdlm("X.txt")

# include("scTenifoldNet.jl")
include("E:\\GitHub\\scTenifoldNet.jl\\src\\scTenifoldNet.jl")
using .scTenifoldNet
A=scTenifoldNet.pcnet(X',3,scalein=true,scaleout=true)

