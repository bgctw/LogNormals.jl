module LogNormals

export plusTwo, AbstractMoments, Moments, QuantilePoint, QuantileSet, 
    @qp, @qp_ll, @qp_l, @qp_m, @qp_u, @qp_uu, 
    @qs_cf90, @qs_cf95


using Distributions, StaticArrays
import DataStructures

"""
    AbstractMoments

A representation of statistical moments of a distribution

The following functions are supported
- length(M::AbstractMoments): get the number of defined moments

The following getters return a single moment or 
throw an error if the moment has not been defined
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
Distributions.mean(M::AbstractMoments) = length(M) >= 1 ? M[1] : 
    error("mean not defined")
Distributions.var(M::AbstractMoments) = length(M) >= 2 ? M[2] : 
    error("variance not defined")
Distributions.std(M::AbstractMoments) = length(M) >= 2 ? sqrt(M[2]) : 
    error("std not defined")
Distributions.skewness(M::AbstractMoments) = length(M) >= 3 ? M[3] : 
    error("skewness not defined")
Distributions.kurtosis(M::AbstractMoments) = length(M) >= 4 ? M[4] : 
    error("kurosis not defined")

struct Moments{N,T} <: AbstractMoments{N}
    all::SVector{N,T}
end
Moments(x...) = Moments(SVector{length(x)}(promote(x...)))
Moments() = Moments(SA[])
Base.getindex(M::Moments, i) = length(M) >= i ? M.all[i] : 
    error("$(i)th moment not defined.")
#Base.convert(::Type{AbstractArray}, M::Moments) = M.all

function moments(D::LogNormal, ::Val{N} = Val(3)) where N 
    N > 4 && error("getting moments above 4 not yet implemented.")
    N == 4 && return(Moments(mean(D), var(D), skewness(D), kurtosis(D)))
    N == 3 && return(Moments(mean(D), var(D), skewness(D)))
    N == 2 && return(Moments(mean(D), var(D)))
    N == 1 && return(Moments(mean(D)))
    N == 0 && return(Moments())
end

function moments(D::Normal, ::Val{N} = Val(3)) where N 
    N > 4 && error("getting moments above 4 not yet implemented.")
    N == 4 && return(Moments(mean(D), var(D), skewness(D), kurtosis(D)))
    N == 3 && return(Moments(mean(D), var(D), skewness(D)))
    N == 2 && return(Moments(mean(D), var(D)))
    N == 1 && return(Moments(mean(D)))
    N == 0 && return(Moments())
end

"""
    fit(Distribution, m::Moments)
    
Fit a statistical distribution to its moments 

# Arguments
- `Distribution`: The type of distribution to fit
- `m`: The moments of the distribution

# Notes
This can be used to approximate one distribution by another

# Examples
```julia
D = fit(LogNormal, Moments(3.2,4.6))
D = fit(LogNormal, moments(Normal(3,1)))
plot(D); lines(!Normal(3,1))
```
"""
function Distributions.fit(::Type{LogNormal}, m::AbstractMoments)
    # https://en.wikipedia.org/wiki/Log-normal_distribution
    length(m) >= 2 || error("Need mean and variance to estimate lognormal")
    γ = 1+var(m)/mean(m)^2
    μ = log(mean(m)/sqrt(γ))
    σ = sqrt(log(γ))
    return LogNormal(μ,σ)
end

struct QuantilePoint
    p
    q
    QuantilePoint(p,q) = 0 < p < 1 ? new(p,q) : 
        error("p must be in (0,1)")
end
QuantilePoint(qp::QuantilePoint; p = qp.p, q = qp.q) = QuantilePoint(p,q)
Base.show(io::IO, qp::QuantilePoint) = print(io, "(p=$(qp.p),q=$(qp.q))")
function Base.isless(x::QuantilePoint,y::QuantilePoint)
    is_equal_q = (x.q == y.q)
    ((x.p == y.p) && !is_equal_q) && error("incompatible QuantilePoints: $x,$y")
    isless = (x.q < y.q)
    # for different p, q needs to be different
    (isless && (x.p > y.p)) && error("incompatible QuantilePoints: $(x),$(y)")
    (!isless && !is_equal_q && (x.p < y.p))  && error("incompatible QuantilePoints: $x,$y")
    return(isless)
end

QuantileSet = DataStructures.SortedSet{QuantilePoint} 

macro qp(p,q) :(QuantilePoint($p, $q)) end
macro qp_ll(q0_025) :(QuantilePoint(0.025, $q0_025)) end
macro qp_l(q0_05) :(QuantilePoint(0.05, $q0_05)) end
macro qp_m(median) :(QuantilePoint(0.5, $median)) end
macro qp_u(q0_95) :(QuantilePoint(0.95, $q0_95)) end
macro qp_uu(q0_975) :(QuantilePoint(0.975, $q0_975)) end

macro qs_cf90(q0_05,q0_95) 
    :(Set([QuantilePoint(0.05,$q0_05),QuantilePoint(0.95,$q0_95)])) end
macro qs_cf95(q0_025,q0_975) 
    :(Set([QuantilePoint(0.025,$q0_025),QuantilePoint(0.975,$q0_975)])) end    

"""
fit(Distribution, qset:QuantileSet)

Fit a statistical distribution to a set of quantiles 

# Arguments
- `Distribution`: The type of distribution to fit
- `qset`: A ordered set of QuantilePoints (p,q)

# Notes
Several macros help to construct ordered sets of QuantilePoints
TODO

# Examples
```julia
DN = fit(Normal, @qp_m(3), @qp_uu(5))
D = fit(LogNormal, @qp_m(3), @qp_uu(5))
quantile.(D, [0.5, 0.975]) ≈ [3,5]
```
"""
function Distributions.fit(::Type{LogNormal}, lower::QuantilePoint, upper::QuantilePoint)
    #length(qset) == 2 || error("only implemented yet for exactly two quantiles.")
    #qset_log = [QuantilePoint(qp, q = log(qp.q)) for qp in qset]
    lower_log = QuantilePoint(lower, q = log(lower.q))
    upper_log = QuantilePoint(upper, q = log(upper.q))
    DN = fit(Normal, lower_log, upper_log)
    LogNormal(params(DN)...)
end

function Distributions.fit(::Type{Normal}, lower::QuantilePoint, upper::QuantilePoint)
    # https://www.johndcook.com/quantiles_parameters.pdf
    if (upper < lower) 
        lower,upper = (upper,lower)
    end
    qz1 = quantile(Normal(), lower.p)
    qz2 = quantile(Normal(), upper.p)
    dqz = (qz2 - qz1)
    σ = (upper.q - lower.q)/dqz
    μ = (lower.q*qz2 - upper.q*qz1)/dqz
    Normal(μ,σ)
end
    
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
