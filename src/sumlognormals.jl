

# function length_itr(x)
#     typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
#         return(length(x))
#     count(x -> true, x)
# end

"""
    sum(dv::AbstractDistributionVector{D}; skipmissings::Val{B} = Val(false))

Compute the distribution of the sum of correlated random variables.

# Arguments
- `dv`: The vector of distributions, see [AbstractDistributionVector](@ref)
- `skipmissing`: Set to Val(true) to conciously care for missing in dv

Optional second arguments are supported 
- correlation matrix, 
- acf::AbstractVector: coefficients of the autocorrelation function starting 
  from lag one

# Examples
```julia
julia> five = plusTwo(3)
5
```
"""
function Base.sum(dv::AbstractDistributionVector{<:D}; 
    skipmissings::Val{B} = Val(false)) where {D,B} end

"""
    sum(dv::AbstractDistributionVector{<:LogNormal})

In addition to `sum(AbstractDistributionVector{<:Distribution})` supports a third 
argument of type `AbstractVector{Bool}` of the same length as `dv`. Flagged
records contribute to the estimate mean of the sum, but not to the decrease
of spread with increasing number of observations.
"""
function Base.sum(dv::DVM; skipmissings::Val{B} = Val(false)) where 
    DVM <: Union{Base.SkipMissing{DV},DV} where 
    {T, DV <: AbstractDistributionVector{LogNormal{T}}, B} 
    B == true && return(sum(skipmissing(dv)))
    # uncorrelated, only sum diagonal
    Ssum = s = zero(T)
    nterm = 0
    for d in dv
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

# own method with argument isgapfilled, because now cannot use
# skipmissing any more and need to allocate nonmissing to 
function Base.sum(dv::DV, isgapfilled::AbstractVector{Bool}; 
    skipmissings::Val{B} = Val(false)) where 
    {T<:Real, DV <: AbstractDistributionVector{LogNormal{T}}, B} 
    length(dv) == length(isgapfilled) || error(
        "argument gapfilled must have the same length as dv ($(length(dv))" *
        "but was $(length(isgapfilled)).")
    if B == true
        # need to allocate anyway with subsetting
        nonmissing = findall(.!ismissing.(dv))
        return(sum(
            @inbounds(dv[nonmissing]), @inbounds(isgapfilled[nonmissing])))
    end
    # uncorrelated, only sum diagonal
    Ssum = s = Ssumnonfilled = zero(T)
    nterm = 0
    for (i,d) in enumerate(dv)
        μ,σ = params(d)
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        if !isgapfilled[i]
            Ssumnonfilled += Si
            s += abs2(σ) * abs2(Si)
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    σ2eff = s/abs2(Ssumnonfilled)
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

function Base.sum(dv::AbstractDistributionVector{<:LogNormal}, 
    acf::AbstractVector{<:Number}, 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)); 
    skipmissings::Val{B} = Val(false), method::Val{S} = Val(:vector)) where 
    {B, S} 
    storage = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(dv))
    if method == Val(:vector) 
        return(sum_lognormals!(
            storage, dv, acf, isgapfilled, skipmissings = skipmissings))
    end
    if method == Val(:bandedmatrix)
        corrM = cormatrix_for_acf(length(dv), acf)
        return(sum_lognormals!(
            storage, dv, corrM, isgapfilled,skipmissings = skipmissings))
    end
    error("Unknown method $method")
end

function sum_lognormals!(S::Vector{Union{Missing,T}}, dv::DV, 
    acf::AbstractVector, 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)); 
    skipmissings::Val{B} = Val(false)) where 
    {DV <: AbstractDistributionVector{D}, B} where D<:LogNormal where T
    #details<< Implements estimation according to
    # Messica A(2016) A simple low-computation-intensity model for approximating
    # the distribution function of a sum of non-identical lognormals for
    # financial applications. 10.1063/1.4964963
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    corrlength = length(acf)
    acfm = vcat(reverse(acf), 1, acf)
    n = length(μ)
    @. S = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing.(S))
    skipmissings != Val(true) && nmissing != 0 && error(
        "Found missing values. Use argument 'skipmissings = Val(true)' " *
        "to sum over nonmissing.")
    # 0 in S has the effect of not contributing to Ssum nor s 
    # For excluding terms of gapfilled records, must make 
    # to either set S[gapfilled] to zero or check in the product
    Spure = disallowmissing(replace(S, missing => 0.0))
    σpure = disallowmissing(replace(σ, missing => 0.0))
    Ssum = sum(Spure)
    # after computing Ssum, can also set gapfilled to zero
    Spure[isgapfilled] .= 0.0
    s = Ssumunfilled = zero(eltype(σpure))
    for i in 1:n
        iszero(Spure[i]) && continue # nothing added
        Ssumunfilled += Spure[i]
        jstart = max(1, i - corrlength)
        jend = min(n, i + corrlength)
        for j in jstart:jend
            acf_ind = (j-i + corrlength +1)
            # sij will be zero if sigma or S is missing (replaced by zero)
            # Sj moved to start to help multiplication by early zero
            sij = Spure[j] * acfm[acf_ind] * σpure[i] * σpure[j] * Spure[i] 
            s += sij
        end
    end
    σ2eff = s/abs2(Ssumunfilled)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, n - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

function Base.sum(dv::AbstractDistributionVector{D}, 
    corr::AbstractMatrix, 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)); 
    skipmissings::Val{B} = Val(false)) where {B, D<:LogNormal}
    S = Vector{Union{Missing,eltype(D)}}(undef, length(dv))
    sum_lognormals!(S, dv, corr, isgapfilled, skipmissings = skipmissings)
end

function sum_lognormals!(S, dv, corr::AbstractMatrix, 
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)); 
    skipmissings::Val{l} = Val(false)) where l
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    # S = allowmissing(similar(μ))
    @. S = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing, S)
    anymissing = nmissing != 0
    l != true && anymissing && error(
         "Found missing values. Use argument 'skipmissings = Val(true)' " *
         "to sum over nonmissing.")
    Ssum::nonmissingtype(eltype(S)) = sum(skipmissing(S))
    # gapfilled records only used for Ssum, can set the to 0 now
    # so they do not contribute to s and Ssumfilled for computation of σ2eff
    S[isgapfilled] .= 0
    Ssumunfilled::nonmissingtype(eltype(S)) = sum(skipmissing(S))
    @. S = σ * S  # do only after Ssum
    # setting S to zero results in summing zero for missing records
    # which is the same as filtering both S and corr
    anymissing && replace!(S, missing => 0.0)
    #s = transpose(disallowmissing(S)) * corr * disallowmissing(S)
    #Spure = view_nonmissing(S) # non-allocating
    Spure = disallowmissing(S) # allocating - tested: is faster than the view
    s = transpose(Spure) * corr * Spure
    σ2eff = s/abs2(Ssumunfilled)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, length(S) - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

