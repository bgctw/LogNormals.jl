"""
    AutoCorrelationFunction{T}

A representation of the autocorrelation function.

It supports accessing the coeficients starting from lag 1 by
- `coef(acf::AutoCorrelationFunction)`: implements StatsBase.coef

Wrapping the vector of coefficients into its own type helps avoiding
method ambiguities.

# Examples
```jldoctest am; output = false, setup = :(using LogNormals)
using StatsBase: coef
acf = AutoCorrelationFunction([0.4,0.1])
coef(acf) == [0.4,0.1]
# output
true
```
"""
struct AutoCorrelationFunction{T}
    coef::T
end
# implements StatsBase coef
coef(acf::AutoCorrelationFunction) = acf.coef
AutoCorrelationFunction(coef::AbstractVector{<:Number}) = 
    AutoCorrelationFunction{typeof(coef)}(coef)

cormatrix_for_acf(n::Int,acf::AutoCorrelationFunction) =
    cormatrix_for_acf(n, coef(acf))

function cormatrix_for_acf(n::Int,acf::AbstractVector) 
    nacf::Int = length(acf)
    corrM = BandedMatrix{Float64}(undef, (n,n), (nacf,nacf))
    corrM[band(0)] .= 1
    for i in 1:nacf
      corrM[band(i)] .= corrM[band(-i)] .= acf[i]
    end
    corrM
end



function sum(dv::AbstractDistributionVector{<:Normal}; 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)), 
    skipmissings::Val{B} = Val(false)) where B
    length(dv) == length(isgapfilled) || error(
        "argument gapfilled must have the same length as dv ($(length(dv))" *
        "but was $(length(isgapfilled)).")
    if B == true
        # need to allocate anyway with subsetting
        nonmissing = findall(.!ismissing.(dv))
        !isempty(nonmissing) && return(sum(
            @inbounds(dv[nonmissing]), 
            isgapfilled = @inbounds(isgapfilled[nonmissing])))
    end
    # uncorrelated, only sum diagonal
    Ssum = s = Ssumnonfilled = zero(eltype(nonmissingtype(eltype(dv))))
    nterm = 0
    for (i,d) in enumerate(dv)
        μ,σ = params(d)
        Ssum += μ
        if !isgapfilled[i]
            Ssumnonfilled += μ
            s += abs2(σ)
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    relerr = √s/Ssumnonfilled
    Normal(Ssum, Ssum * relerr)
end

function sum(dv::AbstractDistributionVector{D}, 
    acf::AutoCorrelationFunction;
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal, DS<:eltype(D), SM}
    # currently use Matrix-method, Maybe implement fast with loop
    corrM = Symmetric(cormatrix_for_acf(length(dv), coef(acf)))
    return(sum_normals(
        dv, corrM, 
        isgapfilled=isgapfilled,skipmissings = skipmissings, storage = storage))
end


function sum(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix}; 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal, DS<:eltype(D), SM}
    sum_normals(
        dv, corr, 
        isgapfilled=isgapfilled, storage=storage, skipmissings=skipmissings)
end

function sum_normals(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix}; 
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
        Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal,SM, DS<:eltype(D)}
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    # var_sum (s) is the sum across all Sigma, i.e. σT * corr * σ
    # missings and gapfilled values do not count -> set to zero
    if SM == true
        #@. storage = (x -> ismissing(x) ? zero(DS) : x)(σ)  
        @. storage = coalesce(σ, zero(DS))  
        Spure = disallowmissing(storage) 
    else
        Spure = disallowmissing(σ)
    end 
    Spure[isgapfilled] .= zero(DS)
    s = transpose(Spure) * corr * Spure
    Ssum = Ssumnonfilled = zero(DS)
    nterm = 0
    for (i,d) in enumerate(dv)
        ismissing(d) && continue
        μ,σ = params(d)
        Ssum += μ
        if !isgapfilled[i]
            Ssumnonfilled += μ
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    relerr = √s/Ssumnonfilled
    #@show Ssum, Ssumnonfilled, s, relerr
    Normal(Ssum, Ssum * relerr)
end



