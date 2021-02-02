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

D = fit(LogNormal; mean = 4, mult_std = exp(1))
exp(params(D)[2])

@test mean(D) ≈ 4 && mult_std(D) ≈ exp(1)


