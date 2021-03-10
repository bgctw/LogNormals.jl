# Strategies for dealing with missings 

The following types support methods to deal with presence of missing values
in their inputs.

```@docs
MissingStrategy
```

Usually, methods should return `missing` if there are
missings present in the input. This `PassMissing()` strategy ensures, 
that missing are not accidently bypassed unconsciously. 
Hence, `PassMissing()` is a good default value for the strategy.

If a caller knows that there might by missing values present, he can consciously
request to handle missings by supplying one of subtypes of `HandleMissingStrategy`.
For an explicit dealing with missings that may cause some cost in computation
use `ExactMissing`. For ignoring missing values use `SkipMissing`.

## Example: Autocorrelation

When dealing with autocorrelation, the position of records in the series
is meaningful. Computations on a series where missing have just been skipped are biased 
because distant records may have become neighbors.

The following example computes the effective number of observations, which tells
how many independent observations would lead to the same uncertainty of the mean. 
It assumes that the errors ``\epsilon_i`` from 
 observation ``x_i = \mu_i + \epsilon_i`` are generated
 by the same process and have the same uncertainty.

```julia
# acf: coefficients of the autocorrelation functions starting from lag 0 
function effective_n_cor(x, acf::AbstractVector, ms::MissingStrategy=PassMissing())
    ms === PassMissing() && Missing <: eltype(x) && any(ismissing.(x)) && return(missing)
    n = length(x)
    k = Base.OneTo(min(n,length(acf))-1) # acf starts with lag 0
    if ms === ExactMissing() && (Missing <: eltype(x))
        # count the number of pairs with missings for each lag
        mk = count_forlags((x_i,x_iplusk)->ismissing(x_i) || ismissing(x_iplusk), x, k)
        nf = n - count(ismissing.(x))
        neff = nf/(1 + 2/nf*sum((n .- k .-mk) .* acf[k.+1]))  
    else
        neff = n/(1 + 2/n*sum((n .- k) .* acf[k.+1]))  
    end
end
```

The default MissingStrategy is `PassMissing()` ensuring that the caller has to
think about how she wants to deal with missings to get proper results.

The formula for the standard error of the mean sums over all covariances, and 
hence the formula can account for missings by counting how many covariance
terms are missing. This counting scans the series for each lag, i.e distance, between
records and therefore comes at a cost.

The caller can decide whether to use the costly exact number of observations by
passing `ExactMissing()` or 
to use a slighly positively biased estimate computed at a lower cost by passing
`SkipMissing()`, for instance when she knows that the potential number of 
missings is small.

## Type stability
Checking the type of inputs for allowing `Missing` in addition to checking for occurence
of missings helps to generate type-stable code.

In the example above: Including the check on `Missing <: eltype(x)`  
leads to a inferred return type of `Float64`
for `x::AbstractVector{Float64}`, 
whereas omitting this check would yield `Union{Missing,Float64}`.

## Design choices
The design strategies alternatively could be passed by keyword arguments,
e.g. `skipmissing::Bool=true` or by an enum type.
Modelling the strategies by a type system allows the dispatch system and the 
compiler to work on them. 

This comes with several advantages:
- Avoid method matching ambiguities across existing methods that cannot
  handle missing values and their extensions by the same name and argument types aside
  from the MissingStrategy.
- Generate efficient code, because `if` branches can be pruned at compile time
- Extensible to further strategies
  
  
