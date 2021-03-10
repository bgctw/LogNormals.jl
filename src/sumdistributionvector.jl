"""
    sum(dv::AbstractDistributionVector, ms=PassMissing(); <keyword arguments>)
    sum(dv::AbstractDistributionVector, corr, ms=PassMissing(); <keyword arguments>)
    sum(dv::AbstractDistributionVector, acf, ms=PassMissing(); <keyword arguments>)
    
    mean(dv::AbstractDistributionVector, ms=PassMissing(); <keyword arguments>)
    mean(dv::AbstractDistributionVector, corr, ms=PassMissing(); <keyword arguments>)
    mean(dv::AbstractDistributionVector, acf, ms=PassMissing(); <keyword arguments>)

Compute the distribution of the sum or mean of correlated random variables.

# Arguments
- `dv`: The vector of distributions, see [`AbstractDistributionVector`](@ref)
- `ms`: [`MissingStrategy`]: set to `SkipMissing()` to consciously care
  for missing values in `dv`.

An optional second arguments supports correlation between random variables.
- `corr::Symmetric`: correlation matrix, or
- `acf::AutoCorrelationFunction`: coefficients of the 
   [`AutoCorrelationFunction`](@ref)

Keyword arguments:
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
end, 
function mean(dv::AbstractDistributionVector)
    error("mean not defined yet for " * 
    "Distributionvector{$(nonmissingtype(eltype(dv)))}")
end



"""
    AutoCorrelationFunction{T}

A representation of the autocorrelation function.

It supports accessing the coeficients starting from lag 0 by
- `coef(acf::AutoCorrelationFunction)`: implements StatsBase.coef
- `coef(acf::AutoCorrelationFunction, lag::Integer)`: coefficient for lag

Wrapping the vector of coefficients into its own type helps avoiding
method ambiguities.

# Examples
```jldoctest am; output = false, setup = :(using LogNormals)
using StatsBase: coef
acf = AutoCorrelationFunction([1,0.4,0.1])
coef(acf) == [1,0.4,0.1]
coef(acf,1) == 0.4
# output
true
```
"""
struct AutoCorrelationFunction{T}
    coef::T
end
# implements StatsBase coef
AutoCorrelationFunction(coef::AbstractVector{<:Number}) = 
    AutoCorrelationFunction{typeof(coef)}(coef)
coef(acf::AutoCorrelationFunction) = acf.coef
coef(acf::AutoCorrelationFunction, lag::Integer) = acf.coef[lag+1]

cormatrix_for_acf(n::Int,acf::AutoCorrelationFunction) =
    cormatrix_for_acf(n, coef(acf))

function cormatrix_for_acf(n::Int,acf::AbstractVector) 
    nacf::Int = length(acf)
    corrM = BandedMatrix{Float64}(undef, (n,n), (nacf-1,nacf-1))
    corrM[band(0)] .= acf[1]
    for i in 1:(nacf-1)
      corrM[band(i)] .= corrM[band(-i)] .= acf[i+1]
    end
    corrM
end
