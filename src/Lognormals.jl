module Lognormals

export plusTwo, AbstractMoments, Moments
using Distributions, StaticArrays

"""
    AbstractMoments

A representation of statistical moments of a distribution

The following functions are supported
- length(M::AbstractMoments): get the number of defined moments

The following getters return a single moment or missing if not defined
- mean(M::AbstractMoments): get the mean
- var(M::AbstractMoments): get the variance
- skewness(M::AbstractMoments): get the variance
- kurtosis(M::AbstractMoments): get the variance
- getindex(M::AbstractMoments,i): get the ith moment, i.e. indexin [i]

Moments-objects be converted to AbstractArray by
- convert(AbstractArray, M::AbstractMoments)

The basic implementation 'Moments' is immutable and
convert(AbstractArray, M::Moments) returns an SArray{N,T}

# Examples
```julia
M = Moments(1,0.2)
length(M) == 2
var(M) == M[2]
ismissing(kurtosis(M))
typeof(convert(AbstractArray,M)) <: AbstractArray
```
"""
abstract type AbstractMoments{N} end
Base.length(::Type{<:AbstractMoments{N}}) where N = N
Base.length(M::AbstractMoments{N}) where N = N
Distributions.mean(M::AbstractMoments) = length(M) >= 1 ? M[1] : missing
Distributions.var(M::AbstractMoments) = length(M) >= 2 ? M[2] : missing
Distributions.skewness(M::AbstractMoments) = length(M) >= 3 ? M[3] : missing
Distributions.kurtosis(M::AbstractMoments) = length(M) >= 4 ? M[4] : missing

struct Moments{N,T} <: AbstractMoments{N}
    all::SVector{N,T}
end
Moments(x...) = Moments(SVector{length(x)}(promote(x...)))
Moments() = Moments(SA[])
Base.getindex(M::Moments, i) = length(M) >= i ? M.all[i] : missing
Base.convert(::Type{AbstractArray}, M::Moments) = M.all

"""
    plusTwo(x)

Sum the numeric "2" to whatever it receives as input

A more detailed explanation can go here, although I guess it is not needed in this case

# Arguments
* `x`: The amount to which we want to add 2

# Notes
* Notes can go here

# Examples
```julia
julia> five = plusTwo(3)
5
```
"""
plusTwo(x) = return x+2

end # module
