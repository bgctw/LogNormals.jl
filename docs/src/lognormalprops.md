# Properties of the LogNormal distribution

The LogNormal distribution can be characterized by
the exponent of its parameters:

- exp(μ): the median
- exp(σ): the multiplicative standard deviation ``\sigma^*``.

Function [`σstar`](@ref) returns the multiplicative standard deviation.

A distribution can be specified by taking the log of median and ``\sigma^*``

```jldoctest; output = false, setup = :(using Distributions, LogNormals)
d = LogNormal(log(2), log(1.2))
σstar(d) == 1.2
# output
true
```

Alternatively the distribution can be specified by its mean and ``\sigma^*`` using type [`Σstar`](@ref)

```jldoctest; output = false, setup = :(using Distributions, LogNormals)
d = fit(LogNormal, 2, Σstar(1.2))
(mean(d), σstar(d)) == (2, 1.2)
# output
true
```

## Detailed API

```@docs
LogNormals.σstar(::LogNormal)
```

```@docs
Distributions.fit(::Type{LogNormal}, ::Any, ::Σstar) 
```

```@docs
LogNormals.Σstar
```
