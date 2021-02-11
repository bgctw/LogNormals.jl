

function length_itr(x)
    typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
        return(length(x))
    count(x -> true, x)
end

# document and example skipmissing(ds)
function sum_lognormals(ds) 
    nterm = length_itr(ds)
    # uncorrelated, only sum diagonal
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    Ssum = s = zero(eltype(eltype(ds)))
    for i in 1:nterm
        μ,σ = params(ds[i])
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        s += abs2(σ) * abs2(Si)
    end
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end

nparams(::Type{<:LogNormal}) = 2
nparams(::Type{<:Normal}) = 2
nparams(::Type{<:LogitNormal}) = 2

abstract type AbstractDistributionSequence{D <: Distribution} end
Base.eltype(::Type{<:AbstractDistributionSequence{D}}) where D = D
function Base.iterate(ds::DS) where DS <: AbstractDistributionSequence
    (ds[1], 1)
end
function Base.iterate(ds::DS, i) where DS <: AbstractDistributionSequence
    i >= length(ds) ? nothing : (ds[i+1], i+1)
end
    
struct DistributionSequence{D <: Distribution, T, Npar} <: AbstractDistributionSequence{D} 
    #params::AbstractMatrix{Union{Missing,T}}
    params::AbstractMatrix
end
DistributionSequence(::Type{D}, parameters) where D <: Distribution =
       DistributionSequence{D, eltype(D), nparams(D)}(parameters)
function DistributionSequence(::Type{D}, pvec::Vararg{Union{Missing,eltype(D)}}) where {D<:Distribution} 
    parms = transpose(hcat(pvec...))
    DistributionSequence(D, parms)
end
function DistributionSequence(ds::Vararg{Union{Missing,D},N}) where {D <: Distribution, N} 
    N == 0 && error(
        "Provide at least one distribution in DistributionSequence(x...).")
    #parms = vcat((transpose.(collect.(params.(ds))))...)
    params_collected = passmissing(collect).(passmissing(params).(ds))
    params_missing = fill(missing, nparams(D))
    parms = hcat((coalesce.(params_collected, Ref(params_missing)))...)
    DistributionSequence(D, parms)        
end

Base.length(ds::DistributionSequence)  =  size(ds.params,2)::Int
function Base.getindex(ds::DistributionSequence{D,T,Npar},i::Int) where {D,T,Npar}
    params_i::SVector{Npar,Union{Missing,T}} = SVector{Npar,Union{Missing,T}}(@view ds.params[:,i]) 
    #params_i::Vector{Union{Missing,T}} = @view ds.params[:,i]
    any(ismissing.(params_i)) && return(missing)
    D(params_i...)
end
StatsBase.params(ds::DS) where DS <: DistributionSequence = ds.params


function Base.sum(ds::DS; skipmissings::Val{B} = Val(false)) where 
    {DS <: Union{AbstractDistributionSequence{LogNormal{T}},
    Base.SkipMissing{AbstractDistributionSequence{LogNormal{T}}}}, B} where T
    if skipmissings == Val(true) 
        return(sum(skipmissing(ds)))
    end
    nterm = length(ds)
    # uncorrelated, only sum diagonal
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    Ssum = s = zero(eltype(eltype(ds)))
    for d in ds
        μ,σ = params(d)
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        s += abs2(σ) * abs2(Si)
    end
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end