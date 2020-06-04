using DelimitedFiles
cd(dirname(@__FILE__))

X=readdlm("X.txt")
Y=readdlm("Y.txt")

# include("scTenifoldNet.jl")
include("E:\\GitHub\\scTenifoldNet.jl\\src\\scTenifoldNet.jl")
using .scTenifoldNet
@time d=scTenifoldNet.manialn(X,Y)
# FC,p,q=scTenifoldNet.drgenes(d)