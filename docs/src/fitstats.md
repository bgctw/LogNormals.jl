# Distribution Fitting to aggregate statistics

This package provides method to fit a distribution to a given
set of aggregate statistics.

```julia
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
```

## Fit to statistical moments

```@docs
Distributions.fit(::Type{<:Distribution}, m::AbstractMoments)
```

```@docs
moments(d::Distribution, ::Val{N} = Val(2)) where N 
```

The syntax `Moments(mean,var)` produces an object of type `Moments <: AbstractMoments`.

```@docs
 AbstractMoments{N}
```

## Fit to several quantile points

```@docs
Distributions.fit(::Type{<:Distribution}, lower::QuantilePoint, upper::QuantilePoint)
```

## Fit to mean,mode,median and a quantile point

```@docs
Distributions.fit(D::Type{<:Distribution}, val, qp::QuantilePoint, ::Val{stats}) where stats
```


