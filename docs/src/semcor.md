# Standard error of the mean of correlated series

The standard error of the mean of a series with positive autocorrelation
is larger than 
for an uncorrelated series, because the variance 
[includes positive covariance
terms](https://en.wikipedia.org/wiki/Variance#Sum_of_correlated_variables) 
in addition to the variance of each record.

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

There are only ``n_{eff}`` number of effective obserservations
in the series.
```@docs
effective_n_cor(x, acf::AbstractVector, ms::MissingStrategy)
```

## Standard error of the mean of a correlated series
```@docs
sem_cor(x::Any, acfe::Any, ms::MissingStrategy; neff=nothing)
```

The default estiamtes the empirical autocorrelation from the given series. 
If possible, use a more precise estimate from longer series. For example
when computing the daily means of an hourly time series, estimate the 
empirical autocorrelation from monthly or annual series and provide it to
the daily applications of `sem_cor` using argument `acfe`.

## Variance of a correlated series
```@docs
var_cor(x::Any, acfe::Any, ms::MissingStrategy; neff=nothing)
```

## Effective autocorrelation function
```@docs
autocor_effective(x, acf::AbstractVector)
```

## Autocorrelation of a series with missing values
```@docs
autocor(x::AbstractVector{Union{Missing,T}}, 
    ms::MissingStrategy=PassMissing(); kwargs... ) where {T<:Real}
```

```@docs
count_forlags(pred::Any, x::Any,lags::AbstractVector) 
```


