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



M = Moments(3)

a = [1,2]

using LogNormals
using StaticArrays, Distributions


D = fit(LogNormal, @qp_m(3), @qp_uu(9))

function Distributions.fit(::Type{LogNormal}; mean::Real, mult_std::Real)
    sigma = log(mult_std)
    mu = log(mean) - sigma*sigma/2
    LogNormal(mu, sigma)
end

mult_std(D::LogNormal) = exp(params(D)[2])

@test mean(D) ≈ 4 && mult_std(D) ≈ exp(1)


function fitstats_docu()
    # to specified moments
    d = fit(LogNormal, Moments(3.0,4.0))
    (mean(d), var(d)) .≈ (3.0, 4.0)

    # to mean and upper quantile point
    d = fit(LogNormal, 3, @qp_uu(8))
    (mean(d), quantile(d, 0.975)) .≈ (3.0, 8.0)
    
    # to mode and upper quantile point
    d = fit(LogNormal, 3, @qp_uu(8), Val(:mode))
    (mode(d), quantile(d, 0.975)) .≈ (3.0, 8.0)

    # to two quantiles, i.e confidence range
    d = fit(LogNormal, @qp_ll(1.0), @qp_uu(8))
    (quantile(d, 0.025), quantile(d, 0.975)) .≈ (1.0, 8.0)

    # approximate a different distribution by matching moments
    dn = Normal(3,2)
    d = fit(LogNormal, moments(dn))
    (mean(d), var(d)) .≈ (3.0, 4.0)
end

# some occurence of true
# https://stackoverflow.com/questions/406230/regular-expression-to-match-a-line-that-doesnt-contain-a-word
occursin.(r"^((?!true|\(true).)*$",["true", "(true, true)", "not true", "other", ""])
# true or (true at the beginning
occursin.(r"^(?!true|\(true).",["true", "(true, true)", "not true", "other", ""])


function experimentAR1()
    using StatsBase
    # generate autocorrelated time series
    using ARFIMA
    N, σ = 1_000, 0.5
    x = arfima(N, σ, nothing, SVector(0.6))
    # re-estimate autocorrelation
    autocor(x, 1:10)
end

tvec = allowmissing([(rand(),rand()) for i=1:6]);
tvec[1] = missing
c = mappedarray(x-> ismissing(x) ? missing : first(x),tvec);
c[1]
c[1:2]
c[[1,2]]


