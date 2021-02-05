# LogNormals.jl

# LogNormals Package

```@docs
 LogNormals
```

see the [github repository](https://github.com/bgctw/LogNormals.jl).


```@meta
DocTestSetup = :(using Pkg; Pkg.add("Distributions"); using Statistics,Distributions,LogNormals)
```

```jldoctest
m = Moments(1,0.2)
n_moments(m) == 2, var(m) == m[2]
# output
(true, true)
```
