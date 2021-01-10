module ScTenifoldNet

using Statistics, LinearAlgebra, Distributions, MultipleTesting, Random, SparseArrays
import TSVD
import TensorToolbox
# import KrylovKit

export pcnet, tenrnet, manialn, drgenes, tensordecomp, sctenifoldnet

const NCOMP1,NCOMP2=3,5
const NLAYERS,NCELLS=10,500

# vecnorm(x) = x./norm.(x[:,i] for i in 1:size(x,2))'
vecnorm(x::AbstractMatrix) = norm.(x[:,i] for i in 1:size(x,2))
function normc!(x)
    for i in 1:size(x,2)
        x[:,i]=x[:,i]./norm(x[:,i])
    end
end

function pcnet(X::AbstractMatrix{T}, p::Int=3;
              scalein::Bool=true, scaleout::Bool=false,
              symmout::Bool=false) where T<:Real
    if scalein
        σ=std(X,dims=1)
        σ[σ.==0].=1.0
        X=(X.-mean(X,dims=1))./σ
    end
    ℊ=size(X,2)
    A=1.0 .-Matrix(I,ℊ,ℊ)
    Threads.@threads for k in 1:ℊ
        y=X[:,k]
        𝒳=X[:,1:end.≠k]
        ϕ=TSVD.tsvd(𝒳,p)[3]
        s=𝒳*ϕ
        s ./= (vecnorm(s).^2)'
        b=sum(y.*s,dims=1)
        𝒷=ϕ*b'
        @inbounds A[k,A[k,:].==1.0]=𝒷
    end
    if symmout
        A=0.5*(A+A')
    end
    if scaleout
        A=A./maximum(abs.(A))
    end
    return convert(Matrix{Float16},A)
  end

function tensordecomp(Λ::AbstractArray{T,3}, p::Int=5;
            scaleout::Bool=true) where T
    𝒯=TensorToolbox.cp_als(Λ,p)
    𝕏=TensorToolbox.full(𝒯)
    A=mean(𝕏[:,:,i] for i=1:size(𝕏,3))
    if scaleout
        A ./=maximum(abs.(A))
    end
    return A
end

function manialn(X::AbstractMatrix{T},Y::AbstractMatrix{T}) where T<:Real
    μ,dim=0.9,30
    n1,n2=size(X,1),size(Y,1)
    W₁,W₂=X.+1,Y.+1
    ℐ=Matrix(I,n1,n2)
    μ = μ*(sum(W₁)+sum(W₂))/(2*sum(ℐ))
    𝕎 = [W₁ μ*ℐ; μ*ℐ' W₂]
    L=diagm(vec(sum(abs.(𝕎),dims=1))).-𝕎
    # λ,V =KrylovKit.eigsolve(L,35,:SR,krylovdim=40)
    # V=hcat(V)
    λ,V = eigen(L)
    i=real(λ).>=1e-8
    V=real(V[:,i])
    dim=min(dim,size(V,2))
    V=V[:,1:dim]
    aln0=V[1:n1,:]
    aln1=V[n1+1:end,:]
    d = norm.((aln0.-aln1)[i,:] for i = 1:n1)
    # _,idx=findmax(dd)
    return d
end

function drgenes(d::AbstractVector{T}) where T<:Real
    d²=d.^2
    fc=d²./mean(d²)
    χ² = Chisq(1)
    pVals = ccdf.(χ², fc)
    pAdjusted = MultipleTesting.adjust(pVals, BenjaminiHochberg())
    return fc,pVals,pAdjusted
end

function tenrnet(X::AbstractMatrix{T}; donorm::Bool=true) where T<:Real
    ℊ,𝒸=size(X)
    if donorm
        lbsz=sum(X,dims=1)
        # X=(X./lbsz)*median(lbsz)
        X=(X./lbsz)*1e4
    end
    A=zeros(Float16, ℊ, ℊ, NLAYERS)
    for k=1:NLAYERS
        println("network ... $k")
        𝕩=X[:,randperm(𝒸)][:,1:NCELLS]    # jackknife (m-out-of-n)
        # 𝕩=X[:,rand(1:𝒸,NCELLS)];            # bootstrapping (m-out-of-n)
        𝕩ᵀ=transpose(𝕩)
        a=pcnet(𝕩ᵀ,NCOMP1)
        a[abs.(a).<quantile(vec(abs.(a)),0.95)].=0.0
        @inbounds A[:,:,k]=sparse(a)
    end
    Z=tensordecomp(A,NCOMP2)
    return Z
end

function sctenifoldnet(X::AbstractMatrix{T}, Y::AbstractMatrix{T}; donorm::Bool=false) where T<:Real
    Z0=tenrnet(X,donorm=donorm)
    Z1=tenrnet(Y,donorm=donorm)
    Z0=0.5*(Z0+Z0')
    Z1=0.5*(Z1+Z1')
    d=manialn(Z0,Z1)
    fc,p,adjp=drgenes(d)
    return d,fc,p,adjp
end

end # module
