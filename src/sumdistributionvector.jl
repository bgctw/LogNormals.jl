"""
    sum(dv::AbstractDistributionVector; <keyword arguments>)
    sum(dv::AbstractDistributionVector, corr; <keyword arguments>)
    sum(dv::AbstractDistributionVector, acf; <keyword arguments>)

Compute the distribution of the sum of correlated random variables.

# Arguments
- `dv`: The vector of distributions, see [`AbstractDistributionVector`](@ref)

An optional second arguments supports correlation between random variables.
- `corr::Symmetric`: correlation matrix, or
- `acf::AutoCorrelationFunction`: coefficients of the 
   [`AutoCorrelationFunction`](@ref)

Keyword arguments:
- `skipmissings`: Set to `Val(true)` to conciously care for missings in dv. 
- `isgapfilled::AbstractVector{Bool}`: set to true for records that should
   contribute to the sum but not to the decrease of relative uncertainty
   with increasing number of records, e.g. for missing records that have
   been estimated (gapfilled). 

The sums of correlated variables require extra allocation and 
support an additional keyword parameter  
- `storage`: a mutable `AbstractVector{eltype(D)}` of length of `dv` 
  that provides storage space to avoid additional allocations.
"""
function sum(dv::AbstractDistributionVector)
    error("sum not defined yet for " * 
    "Distributionvector{$(nonmissingtype(eltype(dv)))}")
end#,
# function sum(dv::AbstractDistributionVector, acf::AutoCorrelationFunction)
#     error("sum with autocorrelation function not defined yet for " * 
#     "Distributionvector{$(nonmissingtype(eltype(dv)))}")
# end,
# function sum(dv::AbstractDistributionVector{D}, 
#     corr::Symmetric) where
#     {D<:Distribution, DS<:eltype(D)}
#     error("sum not defined yet for " * 
#     "Distributionvector{$(nonmissingtype(eltype(dv)))}")
# end


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
