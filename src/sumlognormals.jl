

# function length_itr(x)
#     typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
#         return(length(x))
#     count(x -> true, x)
# end

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
struct DistributionSequence{D <: Distribution, Npar, M <: AbstractMatrix} <: AbstractDistributionSequence{D} 
    #params::AbstractMatrix{Union{Missing,T}}
    params::M
    # inner constructor checking 
end
function DistributionSequence(::Type{D}, parameters::M; 
    npar::Val{Npar} = Val(nparams(D))) where 
    {D<:Distribution{}, Npar, M<:AbstractMatrix} 
    Missing <: eltype(M) || error(
        "Element type of parameter matrix must support missings." *
        " Can you apply 'allowmissing(<yourmatrix>)' in invocation?")
    DistributionSequence{D, Npar, M}(parameters)
end
# DistributionSequence(::Type{D}, parameters) where D <: Distribution =
#        DistributionSequence{D, eltype(D), nparams(D)}(parameters)
function DistributionSequence(::Type{D}, pvec::Vararg{Union{Missing,eltype(D)},Npar}) where 
    {D<:Distribution, Npar} 
    parms = allowmissing(transpose(hcat(pvec...)))
    DistributionSequence(D, parms; npar=Val(Npar))
end
function DistributionSequence(ds::Vararg{Union{Missing,D},N}; 
    npar::Val{Npar} = Val(nparams(D))) where {D <: Distribution, N, Npar} 
    N == 0 && error(
        "Provide at least one distribution in DistributionSequence(x...).")
    # not very efficient, but probably used for only few Distributions
    # prefer the other constructors for large sequences.
    # params may be tuples of different type and may hold missings
    params_collected = passmissing(collect).(
        map(passmissing(x -> promote(x...)), (passmissing(params).(ds))))
    params_missing = fill(missing, Npar) # inserted for missing
    parms = allowmissing(hcat((coalesce.(params_collected, Ref(params_missing)))...))
    DistributionSequence(D, parms, npar=npar)        
end

Base.length(ds::DistributionSequence)  =  size(ds.params,2)::Int
function Base.getindex(ds::DistributionSequence{D,Npar,M},i::Int)::Union{Missing, D} where {D,Npar,M}
    params_i::SVector{Npar,Union{Missing,eltype(M)}} = SVector{Npar,Union{Missing,eltype(M)}}(@view ds.params[:,i]) 
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

function cormatrix_for_acf(n::Int,acf::AbstractVector) 
    nacf::Int = length(acf)
    corrM = BandedMatrix{Float64}(undef, (n,n), (nacf,nacf))
    corrM[band(0)] .= 1
    for i in 1:nacf
      corrM[band(i)] .= corrM[band(-i)] .= acf[i]
    end
    corrM
end

function Base.sum(ds::DS, acf::AbstractVector; skipmissings::Val{B} = Val(false), method::Val{S} = Val(:vector)) where 
    {DS <: AbstractDistributionSequence{LogNormal}, B, S} 
    storage = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(ds))
    if method == Val(:vector) 
        return(sum_lognormals!(storage, ds, acf, skipmissings = skipmissings))
    end
    if method == Val(:bandedmatrix)
        corrM = cormatrix_for_acf(length(ds), acf)
        return(sum_lognormals!(storage, ds, corrM, skipmissings = skipmissings))
    end
    error("Unknown method $method")
end


function sum_lognormals!(S::Vector{Union{Missing,T}}, ds::DS, acf::AbstractVector; 
    skipmissings::Val{B} = Val(false)) where 
    {DS <: AbstractDistributionSequence{D}, B} where D<:LogNormal where T
    ##details<< Implements estimation according to
    ## Messica A(2016) A simple low-computation-intensity model for approximating
    ## the distribution function of a sum of non-identical lognormals for
    ## financial applications. 10.1063/1.4964963
    parms = params(ds)
    μ = @view parms[1,:]
    σ = @view parms[2,:]
    corrlength = length(acf)
    acfm = vcat(reverse(acf), 1, acf)
    n = size(parms,2)
    @. S = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing.(S))
    S2 = disallowmissing(S)
    S2[1] = 3.0
    replace!(S, missing => 0.0)
    Ssum = s = 0.0 #zero(eltype(LogNormal))
    for i in 1:n
        Ssum += S[i]
        jstart = max(1, i - corrlength)
        jend = min(n, i + corrlength)
        for j in jstart:jend
            acf_ind = (j-i + corrlength +1)
            sij = acfm[acf_ind] * σ[i] * σ[j] * S[i] * S[j]
            if !ismissing(sij) 
                s += sij
            end
        end
    end
    skipmissings != Val(true) && nmissing != 0 && error(
        "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
    #σ2eff::eltype(LogNormal) = s/abs2(Ssum)
    #μ_sum::eltype(LogNormal) = log(Ssum) - σ2eff/2
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, n - nmissing
    LogNormal(μ_sum, √σ2eff)  
end


function Base.sum(ds::DS, corr::AbstractMatrix; skipmissings::Val{B} = Val(false)) where 
    {DS <: AbstractDistributionSequence{LogNormal}, B} 
    S = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(ds))
    sum_lognormals!(S, ds, corr, skipmissings = skipmissings)
end

view_nonmissing(x) = of_eltype(nonmissingtype(eltype(x)),x)

function sum_lognormals!(S, ds, corr::AbstractMatrix; 
    skipmissings::Val{l} = Val(false)) where l
    parms = params(ds)
    μ = @view parms[1,:]
    σ = @view parms[2,:]
    # S = allowmissing(similar(μ))
    @. S = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing, S)
    anymissing = nmissing != 0
    skipmissings != Val(true) && anymissing && error(
         "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
    Ssum::nonmissingtype(eltype(S)) = sum(skipmissing(S))
    @. S = σ * S  # do only after Ssum
    # setting S to zero results in summing zero for missing records
    # which is the same as filtering both S and corr
    anymissing && replace!(S, missing => 0.0)
    #s = transpose(disallowmissing(S)) * corr * disallowmissing(S)
    #Spure = view_nonmissing(S) # non-allocating
    Spure = disallowmissing(S) # allocating - tested: is faster than the view
    s = transpose(Spure) * corr * Spure
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, length(S) - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

