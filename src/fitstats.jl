"""
    AbstractMoments{N}

A representation of statistical moments of a distribution

The following functions are supported
- `n_moments(m)`: get the number of recorded moments

The following getters return a single moment or 
throw an error if the moment has not been recorded
- `mean(m)`: get the mean
- `var(m)`: get the variance
- `skewness(m)`: get the variance
- `kurtosis(m)`: get the variance
- `getindex(m,i)`: get the ith moment, i.e. indexing m[i]

The basic implementation `Moments` is immutable and
`convert(AbstractArray, m::Moments)` returns an `SArray{N,T}`.

# Examples
```jldoctest am; output = false, setup = :(using Statistics,Distributions,LogNormals)
m = Moments(1,0.2);
n_moments(m) == 2
var(m) == m[2]
# output
true
```
```julia
kurtosis(m) # throws error because its above 2nd moment
```
"""
abstract type AbstractMoments{N} end
n_moments(::Type{<:AbstractMoments{N}}) where N = N
n_moments(m::AbstractMoments{N}) where N = N
Distributions.mean(m::AbstractMoments) = n_moments(m) >= 1 ? m[1] : 
    error("mean not recorded")
Distributions.var(m::AbstractMoments) = n_moments(m) >= 2 ? m[2] : 
    error("variance not recorded")
Distributions.std(m::AbstractMoments) = n_moments(m) >= 2 ? sqrt(m[2]) : 
    error("std not recorded")
Distributions.skewness(m::AbstractMoments) = n_moments(m) >= 3 ? m[3] : 
    error("skewness not recorded")
Distributions.kurtosis(m::AbstractMoments) = n_moments(m) >= 4 ? m[4] : 
    error("kurtosis not recorded")


struct Moments{N,T} <: AbstractMoments{N}
    all::SVector{N,T}
end
Moments(x...) = Moments(SVector{length(x)}(promote(x...)))
Moments() = Moments(SA[])
Base.getindex(m::Moments, i) = n_moments(m) >= i ? m.all[i] : 
    error("$(i)th moment not recorded.")
Base.convert(::Type{AbstractArray}, m::Moments) = m.all

"""
    moments(D, ::Val{N} = Val(2))

Get the first N moments of a distribution.

See also type [`AbstractMoments`](@ref).

## Examples
```jldoctest; output = false, filter = r"Moments{4,", setup = :(using Statistics,Distributions,LogNormals)
moments(LogNormal(), Val(4))  # first four moments 
moments(Normal())  # mean and variance

# output
Moments{2,Float64}([0.0, 1.0])
```
"""
function moments(d::Distribution, ::Val{N} = Val(2)) where N 
    typeof(N) <: Integer || error("N must be a positive Integer")
    N > 4 && error("Getting moments above 4 not yet implemented for distribution $(typeof(d)).")
    N == 4 && return(Moments(mean(d), var(d), skewness(d), kurtosis(d)))
    N == 3 && return(Moments(mean(d), var(d), skewness(d)))
    N == 2 && return(Moments(mean(d), var(d)))
    N == 1 && return(Moments(mean(d)))
    N == 0 && return(Moments())
    error("N must be a positive Integer.")
end

"""
    fit(D, m)
    
Fit a statistical distribution of type `D` to given moments `m`.

# Arguments
- `D`: The type of distribution to fit
- `m`: The moments of the distribution

# Notes
This can be used to approximate one distribution by another.

See also [`AbstractMoments`](@ref), [`moments`](@ref). 


# Examples
```jldoctest fm1; output = false, setup = :(using Statistics,Distributions,LogNormals)
d = fit(LogNormal, Moments(3.2,4.6));
(mean(d), var(d)) .≈ (3.2,4.6)
# output
(true, true)
```
```jldoctest fm1; output = false, setup = :(using Statistics,Distributions,LogNormals)
d = fit(LogNormal, moments(Normal(3,1.2)));
(mean(d), std(d)) .≈ (3,1.2)
# output
(true, true)
```
```julia
plot(d); lines(!Normal(3,1.2))
```
"""
Distributions.fit(::Type{D}, m::AbstractMoments) where {D<:Distribution} = 
    error("fitting to moments not implemented for distribution of type $D")



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

macro qp(p,q) :(QuantilePoint($(esc(p)), $(esc(q)))) end
macro qp_ll(q0_025) :(QuantilePoint(0.025, $(esc(q0_025)))) end
macro qp_l(q0_05) :(QuantilePoint(0.05, $(esc(q0_05)))) end
macro qp_m(median) :(QuantilePoint(0.5, $(esc(median)))) end
macro qp_u(q0_95) :(QuantilePoint(0.95, $(esc(q0_95)))) end
macro qp_uu(q0_975) :(QuantilePoint(0.975, $(esc(q0_975)))) end

macro qs_cf90(q0_05,q0_95) 
    :(Set([QuantilePoint(0.05,$(esc(q0_05))),QuantilePoint(0.95,$(esc(q0_95)))])) end
macro qs_cf95(q0_025,q0_975) 
    :(Set([QuantilePoint(0.025,$(esc(q0_025))),QuantilePoint(0.975,$(esc(q0_975)))])) end    


"""
    fit(D, lower, upper)

Fit a statistical distribution to a set of quantiles 

# Arguments
- `D`: The type of the distribution to fit
- `lower`:  lower QuantilePoint (p,q)
- `upper`:  upper QuantilePoint (p,q)

# Notes
Several macros help to construct QuantilePoints
- `@qp(p,q)`    quantile at specified p: `QuantilePoint(p, q)`
- `@qp_ll(q0_025)`  quantile at very low p: `QuantilePoint(0.025, q0_025)` 
- `@qp_l(q0_05)`    quantile at low p: `QuantilePoint(0.05, q0_05)` 
- `@qp_m(median)`   quantile at median: `QuantilePoint(0.5, median)` 
- `@qp_u(q0_95)`    quantile at high p: `QuantilePoint(0.95, q0_95)`  
- `@qp_uu(q0_975)`  quantile at very high p: `QuantilePoint(0.975, q0_975)` 

# Examples
```jldoctest; output = false, setup = :(using Statistics,Distributions,LogNormals)
d = fit(LogNormal, @qp_m(3), @qp_uu(5));
quantile.(d, [0.5, 0.975]) ≈ [3,5]
# output
true
```
"""
Distributions.fit(::Type{D}, lower::QuantilePoint, upper::QuantilePoint) where D<:Distribution =
    error("fitting to two quantile points not implemented for distribution of type $D")

fit_median_quantile(D::Type{DT}, median, qp::QuantilePoint) where {DT <: Distribution} = 
    return(fit(D, @qp_m(median), qp))



"""
    fit(D, val, qp, ::Val{stats} = Val(:mean))

Fit a statistical distribution to a quantile and given statistics

# Arguments
- `D`: The type of distribution to fit
- `val`: The value of statistics
- `qp`: QuantilePoint(p,q)
- `stats` Which statistics to fit: defaults to `Val(:mean)`. 
   Alternatives are: `Val(:mode)`, `Val(:median)`

# Examples
```jldoctest fm2; output = false, setup = :(using Statistics,Distributions,LogNormals)
d = fit(LogNormal, 5, @qp_uu(14));
(mean(d),quantile(d, 0.975)) .≈ (5,14)
# output
(true, true)
```
```jldoctest fm2; output = false, setup = :(using Statistics,Distributions,LogNormals)
d = fit(LogNormal, 5, @qp_uu(14), Val(:mode));
(mode(d),quantile(d, 0.975)) .≈ (5,14)
# output
(true, true)
```
"""
function Distributions.fit(::Type{D}, val, qp::QuantilePoint, ::Val{stats} = Val(:mean)) where {D<:Distribution, stats}
    stats == :mean && return(fit_mean_quantile(D, val, qp))
    stats == :mode && return(fit_mode_quantile(D, val, qp))
    stats == :median && return(fit_median_quantile(D, val, qp))
    error("unknown stats: $stats")
end;




