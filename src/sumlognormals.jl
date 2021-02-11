

function length_itr(x)
    typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
        return(length(x))
    count(x -> true, x)
end

nparams(::Type{<:LogNormal}) = 2
nparams(::Type{<:Normal}) = 2
nparams(::Type{<:LogitNormal}) = 2

## AbstractDistributionSequence
"""
    AbstractDistributionSequence{D <: Distribution}

Is any type able iterating the same type of Distribution.
Parametrized by a `Distribution` defining the type of the Distribution.

Items may be missing. Hence the element type of the Iterator is `Union{Missing,D}`.

Specific implementations, such as [DistributionSequence](@ref) need
to implement methods `length` and `getindex`.
"""
abstract type AbstractDistributionSequence{D <: Distribution} end
Base.eltype(::Type{<:AbstractDistributionSequence{D}}) where D = Union{Missing,D}
function Base.iterate(ds::DS, state=1) where DS <: AbstractDistributionSequence
    state > length(ds) ? nothing : (ds[state], state+1)
end
Base.iterate(rds::Iterators.Reverse{DS}, state=length(rds.itr)) where 
    DS <: AbstractDistributionSequence = 
    state < 1 ? nothing : (rds.itr[state], state-1)
StatsBase.params(ds::DS) where DS <: AbstractDistributionSequence{D} where D = 
    hcat(collect.(params.(ds))...)


## DistributionSequence   
"""
    DistributionSequence{D <: Distribution, T, Npar}

Is an Matrix-based implementation of [AbstractDistributionSequence](@ref).
It is parametrized by the type of `Distribution` D, the element type of D,
and the number of parameters of D.

Method `params` returns a (Npar x N) AbstractArray of distribution parameters
"""
struct DistributionSequence{D <: Distribution, T, Npar} <: AbstractDistributionSequence{D} 
    #params::AbstractMatrix{Union{Missing,T}}
    params::AbstractMatrix
    # inner constructor checking 
end
DistributionSequence(::Type{D}, parameters; npar::Val{Npar} = Val(nparams(D))) where {D <: Distribution, Npar} =
       DistributionSequence{D, eltype(D), Npar}(parameters)
# DistributionSequence(::Type{D}, parameters) where D <: Distribution =
#        DistributionSequence{D, eltype(D), nparams(D)}(parameters)
function DistributionSequence(::Type{D}, pvec::Vararg{Union{Missing,eltype(D)},Npar}) where {D<:Distribution, Npar} 
    parms = transpose(hcat(pvec...))
    DistributionSequence(D, parms; npar=Val(Npar))
end
function DistributionSequence(ds::Vararg{Union{Missing,D},N}; npar::Val{Npar} = Val(nparams(D))) where {D <: Distribution, N, Npar} 
    N == 0 && error(
        "Provide at least one distribution in DistributionSequence(x...).")
    #parms = vcat((transpose.(collect.(params.(ds))))...)
    params_collected = passmissing(collect).(passmissing(params).(ds))
    params_missing = fill(missing, Npar)
    parms = hcat((coalesce.(params_collected, Ref(params_missing)))...)
    DistributionSequence(D, parms, npar=npar)        
end

Base.length(ds::DistributionSequence)  =  size(ds.params,2)::Int
function Base.getindex(ds::DistributionSequence{D,T,Npar},i::Int)::Union{Missing, D} where {D,T,Npar}
    params_i::SVector{Npar,Union{Missing,T}} = SVector{Npar,Union{Missing,T}}(@view ds.params[:,i]) 
    #params_i::Vector{Union{Missing,T}} = @view ds.params[:,i]
    any(ismissing.(params_i)) && return(missing)
    D(params_i...)
end
StatsBase.params(ds::DS) where DS <: DistributionSequence = ds.params
#StatsBase.params(ds::DS, i::Int) where DS <: DistributionSequence = ds.params[:,i]


function Base.sum(ds::DSM; skipmissings::Val{B} = Val(false)) where 
    DSM <: Union{Base.SkipMissing{DS},DS} where 
    {DS <: AbstractDistributionSequence{LogNormal{T}}, B} where T
    skipmissings == Val(true) && return(sum(skipmissing(ds)))
    # uncorrelated, only sum diagonal
    Ssum = s = zero(T)
    nterm = 0
    for d in ds
        μ,σ = params(d)
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        s += abs2(σ) * abs2(Si)
        nterm += 1
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end

function sum_lognormals!(S, ds, corr::AbstractMatrix; corrlength = length(S)-1, 
    skipmissings::Val{l} = Val(false)) where l
    ##details<< Implements estimation according to
    ## Messica A(2016) A simple low-computation-intensity model for approximating
    ## the distribution function of a sum of non-identical lognormals for
    ## financial applications. 10.1063/1.4964963
    parms = params(ds)
    μ = @view parms[1,:]
    σ = @view parms[2,:]
    Ssum = s = zero(eltype(LogNormal))
    nterm = 0
    n = size(parms,2)
    for i in 1:n
        S[i] = Si = exp(μ[i] + abs2(σ[i])/2)
        if !ismissing(Si)
            nterm += one(nterm)
            Ssum += Si
            jstart = max(1, i - corrlength)
            jend = min(nterm, i + corrlength)
            for j in jstart:jend
                sij = corr[i,j] * σ[i] * σ[j] * S[i] * S[j]
                if !ismissing(sij) 
                    s += sij
                end
            end
        end
    end
    skipmissings != Val(true) && nterm != n && error(
        "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    @show nterm, Ssum, s, σ2eff, μ_sum
    LogNormal(μ_sum, √σ2eff)  
end

function Base.sum(ds::DS, corr::AbstractMatrix; corrlength = size(corr,1)-1, skipmissings::Val{B} = Val(false)) where 
    {DS <: AbstractDistributionSequence{LogNormal}, B} 
    S = Vector{Union{Missing,eltype(LogNormal)}}(undef, size(corr,1))
    sum_lognormals!(S, ds, corr, corrlength = corrlength, skipmissings = skipmissings)
end

