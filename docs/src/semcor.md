# Standard error of the mean of correlated series

## Effective number of observations

For uncorrelated series, the variance of the mean decreases by the number of records compared to the 
variance of a single record.

```math
Var(\bar{x}) = {Var(x) \over n}
```

For correlated series, the effective number of observations is defined as to give the same scaling.

```math
Var(\bar{x}) = {Var(x) \over n_{eff}}
```

Since $n_{eff} < n$ for positive correlations, the 
uncertainty of the mean, i.e. standard deviation of 
the mean is larger than that for an uncorrelated series.

```@docs
effective_n_cor(x, acf::AbstractVector; exactmissing::Bool=true)
```

## Standard error of the mean of a correlated series
```@docs
sem_cor(x::Any, acfe::Any; exactmissing::Bool=true, neff=nothing)
```

## Variance of a correlated series
```@docs
var_cor(x::Any, acfe::Any; exactmissing::Bool=true, neff=nothing)
```

## Effective autocorrelation function
```@docs
autocor_effective(x::Any, acf::Any)
```

## Autocorrelation of a series with missing values
```@docs
 autocor(x::AbstractVector{Union{Missing,T}}, 
    lags::AbstractVector{<:Integer}=StatsBase.default_autolags(size(x,1)); 
    demean::Bool=true, exactmissing::Bool=true, kwargs...) where T <: Real
```

```@docs
count_forlags(pred::Any, x::Any,lags::AbstractVector) 
```


