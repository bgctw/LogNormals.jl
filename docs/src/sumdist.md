# Aggreatation of random variables

## Sum and Mean of several correlated random variables
```@docs
sum(::AbstractDistributionVector)
```

## Helpers
If correlations are only dependent on the distance of records, one can specify
these correlation by a vector starting with distance, i.e. lag, 0.
```@docs
AutoCorrelationFunction
```
