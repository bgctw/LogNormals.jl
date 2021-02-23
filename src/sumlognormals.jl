

# function length_itr(x)
#     typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
#         return(length(x))
#     count(x -> true, x)
# end

"""
    sum(dv::AbstractDistributionVector; skipmissings::Val{B} = Val(false))

Compute the distribution of the sum of correlated random variables.

# Arguments
- `dv`: The vector of distributions, see [`AbstractDistributionVector`](@ref)
- `skipmissing`: Set to `Val(true)` to conciously care for missings in dv. 
   By default missings may result in errors thrown.

Optional second arguments are supported 
- `corr::Symmetric(T, <:AbstractMatrix) where T`: correlation matrix, 
- `acf::AbstractVector`: coefficients of the autocorrelation function starting 
  from lag one

The sums of correlated variables require extra allocation and 
support an additional keyword parameter  
- `storage`: a mutable `AbstractVector{eltype(D)}` of length of `dv` 
  that provides storage space to avoid additional allocations.
"""
function Base.sum(dv::AbstractDistributionVector; 
    skipmissings::Val{B} = Val(false)) where {B} 
    error("sum not defined yet for Distributionvector{$(nonmissingtype(eltype(dv)))}")
end

"""
    sum(dv::AbstractDistributionVector{<:LogNormal})

In addition to [`sum(AbstractDistributionVector{<:Distribution})`](@ref) 
supports a third 
argument of type `AbstractVector{Bool}` of the same length as `dv`. Flagged
records contribute to the estimate mean of the sum, but not to the decrease
of spread with increasing number of observations.
"""
# function Base.sum(dv::DVM; skipmissings::Val{B} = Val(false)) where 
#     {T, DV <: AbstractDistributionVector{LogNormal{T}}, DVM <: Union{Base.SkipMissing{DV},DV}, B} 
# need to keep Type T, to match LogNormal{Float64}, <:LogNormal is also a UnionAll
function Base.sum(dv::Union{Base.SkipMissing{DV},DV}; 
    skipmissings::Val{B} = Val(false)) where 
    {DV <: AbstractDistributionVector{<:LogNormal}, B} 
    B == true && return(sum(skipmissing(dv)))
    # uncorrelated, only sum diagonal
    Ssum = s = zero(eltype(nonmissingtype(eltype(dv))))
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
function Base.sum(dv::AbstractDistributionVector{<:LogNormal}, 
    isgapfilled::AbstractVector{Bool}; 
    skipmissings::Val{B} = Val(false)) where B
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
    Ssum = s = Ssumnonfilled = zero(eltype(nonmissingtype(eltype(dv))))#zero(T)
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

function Base.sum(dv::AbstractDistributionVector{D}, 
    acf::AbstractVector{<:Number}, 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,ST}} = 
        Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false),    
    method::Val{M} = Val(:vector)) where 
    {D<:LogNormal, SM, ST<:eltype(D), M} 
    #storage = Vector{Union{Missing,eltype(D)}}(undef, length(dv))
    if M == :vector
        return(sum_lognormals(
            dv, acf, isgapfilled, storage = storage, skipmissings = skipmissings))
    end
    if M == :bandedmatrix
        corrM = Symmetric(cormatrix_for_acf(length(dv), acf))
        return(sum_lognormals(
            dv, corrM, isgapfilled,skipmissings = skipmissings, storage = storage))
    end
    error("Unknown method $method")
end

function sum_lognormals(dv::AbstractDistributionVector{D}, 
    acf::AbstractVector, 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SK} = Val(false)) where 
    #{D<:LogNormal, ST<:, B}
    {D<:LogNormal, DS<:eltype(D), SK}
    #details<< Implements estimation according to
    # Messica A(2016) A simple low-computation-intensity model for approximating
    # the distribution function of a sum of non-identical lognormals for
    # financial applications. 10.1063/1.4964963
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    corrlength = length(acf)
    acfm = vcat(reverse(acf), 1, acf)
    n = length(μ)
    @. storage = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing.(storage))
    skipmissings != Val(true) && nmissing != 0 && error(
        "Found missing values. Use argument 'skipmissings = Val(true)' " *
        "to sum over nonmissing.")
    # 0 in storage has the effect of not contributing to Ssum nor s 
    # For excluding terms of gapfilled records, must make 
    # to either set storage[gapfilled] to zero or check in the product
    Spure = disallowmissing(replace(storage, missing => 0.0))
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
            # sij will be zero if sigma or storage is missing (replaced by zero)
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
    corr::Symmetric{DS,<:AbstractMatrix}, 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:LogNormal, DS<:eltype(D), SM}
    sum_lognormals(
        dv, corr, isgapfilled, storage = storage, skipmissings = skipmissings)
end

function sum_lognormals(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix}, 
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,DS}} = 
        Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:LogNormal,SM, DS<:eltype(D)}
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    # storage = allowmissing(similar(μ))
    @. storage = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing, storage)
    anymissing = nmissing != 0
    SM != true && anymissing && error(
         "Found missing values. Use argument 'skipmissings = Val(true)' " *
         "to sum over nonmissing.")
    Ssum::nonmissingtype(eltype(storage)) = sum(skipmissing(storage))
    # gapfilled records only used for Ssum, can set the to 0 now
    # so they do not contribute to s and Ssumfilled for computation of σ2eff
    storage[isgapfilled] .= 0
    Ssumunfilled::nonmissingtype(eltype(storage)) = sum(skipmissing(storage))
    @. storage = σ * storage  # do only after Ssum
    # setting storage to zero results in summing zero for missing records
    # which is the same as filtering both storage and corr
    anymissing && replace!(storage, missing => 0.0)
    #s = transpose(disallowmissing(storage)) * corr * disallowmissing(storage)
    #Spure = view_nonmissing(storage) # non-allocating
    Spure = disallowmissing(storage) # allocating - tested: is faster than the view
    s = transpose(Spure) * corr * Spure
    σ2eff = s/abs2(Ssumunfilled)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, length(storage) - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

