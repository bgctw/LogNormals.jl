"""
    count_forlags(pred, x, lags)
    count_forlag(pred, x, k::Integer)

Count the number of pairs for lag `k` which fulfil a predicate.

# Arguments
- `pred::Function(x_i,x_iplusk)::Bool`: The predicate to be applied to each pair 
- `x`: The series whose lags are inspected.
- `lags`: An iterator of Integer lag sizes
- `k`: A single lag.

Common case is to compute the number of missings for the autocorrelation:
with predicate `missinginpair(x,y) = ismissing(x) || ismissing(y)`.

"""
function count_forlag(pred,x,k::Integer)
    s = 0
    ax = OffsetArrays.no_offset_view(axes(x,1))
    for i in 1:(length(x)-k)
        if pred(x[ax[i]], x[ax[i+k]]); s += 1; end
    end
    s
end,
function count_forlags(pred, x,lags::AbstractVector) 
    count_forlag.(tuple(pred), tuple(x),lags)
end

"""
    autocor(x::AbstractVector{x::Union{Missing,<:T}, lags, ms::MS=PassMissing(); dmean::Bool=true}

Estimate the autocorrelation function accounting for missing values.

# Arguments
- `x`: series, which may contain missing values
- `lags`: Integer vector of the lags for which correlation should be computed
- `ms`: [`MissingStrategy`](@ref). If `ExactMissing()` then divide the sum
   in the formula of the exepected value in the formula for the correlation
   at lag `k` by `n - nmissing` instead of `n`, 
   where `nimissing` is the number of records where there is a missing either
   in the original vector or its lagged version (see [`count_forlags`](@ref)).
- `deman`: if `false`, assume `mean(x)==0`.
- `skipmissings`: need to explicitly request `skipmissings=Val(true)` to 
    to consciously deal with missing values.

If the missing strategy is set to `SkipMissing()` then the computation is faster, 
but it is more strongly biased low with increasing number of missings. 
Note that `StatsBase.autocor` uses devision by `n` instead of 'n-k', the
true length of the vectors correlated at lag `k` resulting in 
low-biased correlations of higher lags for numerical stability reasons.
"""
function autocor(x::AbstractVector{Union{Missing,T}}, 
    ms::MissingStrategy=PassMissing(); kwargs... ) where {T<:Real} 
    autocor(x, StatsBase.default_autolags(size(x,1)), ms, kwargs...)
end,
function autocor(x::AbstractVector{Union{Missing,T}}, 
    lags::AbstractVector{<:Integer},
    ms::MS=PassMissing(); demean::Bool=true, kwargs... ) where 
    {T<:Real,MS<:MissingStrategy}
    # if not inheriting from missing, use the StatsBase standard
    !(Missing <: eltype(x)) && return(autocor(x, lags; demean=demean, kwargs...))
    MS == PassMissing && any(ismissing.(x)) && return(missing)
    # only set to zero after demeaning, otherwise correlation estimates increase
    z::Vector{Union{Missing,T}} = demean ? x .- mean(skipmissing(x)) : x
    # replace missing by zero: new type will not match signature of current function 
    zpure = coalesce.(z, zero(z))::Vector{T} 
    acf = autocor(zpure,lags;demean=false,kwargs...)
    MS != ExactMissing && return(acf)
    # correct for sum has been devided by a larger number of terms
    # (including terms of missing/0) 
    lx = length(x)
    #zcorr = nterm_forlag(x,0)/lx
    missinginpair(x,y) = ismissing(x) || ismissing(y)
    zcorr_can = lx - count_forlag(missinginpair,x,0)#nterm_forlag(x,0)
    for (i,k) in enumerate(lags)
        #autocor code devides by lx instead of (lx-k)
        #https://github.com/JuliaStats/StatsBase.jl/issues/273#issuecomment-307560660
        # this cancels with lx of zcorr
        #acf[i] *= (lx - k)/nterm_forlag(x,k) * zcorr
        acf[i] *= zcorr_can/(lx - count_forlag(missinginpair,x,k)) 
    end
    acf
end
# accept MissingStrategy for methods with Missing not part of eltype
autocor(x::StatsBase.RealVector, ms::MissingStrategy; kwargs...) = autocor(x; kwargs...)
autocor(x::StatsBase.RealVector, lags::StatsBase.IntegerVector, ms::MissingStrategy; 
    kwargs...) = autocor(x, lags; kwargs...)

"""
    autocor_effective(x, ms::MissingStrategy=PassMissing())
    autocor_effective(x, acf)

Estimate the effective autocorrelation function for series x.

# Arguments
- `x`: An iterator of a series of observations
- `ms`: [`MissingStrategy`](@ref) passed to [`autocor`](@ref)
- `acf`: AutocorrelationFunction starting from lag 0

# Notes
- The effect autocorrelation function  
  are the first coefficients of the autocorrelation function up to 
  before the first negative coefficient. 
- According to Zieba 2011 using this effective version rather the full version
  when estimating the autocorrelationfunction from the data
  yields better result for the standard error of the mean ([`sem_cor`](@ref)).
- Optional argument `acf` allows the caller to provide a precomputed estimate
  of autocorrelation function (see [`autocor`](@ref)).
"""
function autocor_effective(x, ms::MissingStrategy=PassMissing()) 
    autocor_effective(x, autocor(x, ms))
end,
function autocor_effective(x, acf::AbstractVector)
    #maybe implement a more efficient version that computes only the
    # first lags and further lags if not found negative correlation
    i = findfirst(x -> x <= 0.0, acf)
    isnothing(i) && return(acf)
    acf[1:(i-1)]
end


@doc raw"""
    sem_cor(x, ms::MissingStrategy=PassMissing())
    sem_cor(x, acf::AbstractVector, ms::MissingStrategy=PassMissing())

Estimate the standard error of the mean of an autocorrelated series:
``Var(\bar{x}) = {Var(x) \over n_{eff}}``.    

# Arguments
- `x`: An iterator of a series of observations
- `acf`: AutocorrelationFunction starting from lag 0. 
- `ms`: [`MissingStrategy`](@ref) passed to [`effective_n_cor`](@ref).
  Value of `SkipMissing()` speeds up computation compared to `ExactMissing()`,
  but leads to a negatively biased result with absolute value of the bias 
  increasing with the number of missings.

# Optional Arguments
- `neff`: may provide a precomputed number of observations for efficiency.
"""
function sem_cor(x, ms::MissingStrategy=PassMissing(); kwargs...) 
    sx = typeof(ms) <: HandleMissingStrategy ? std(skipmissing(x)) : std(x)::eltype(x)
    ea = early_var_return(x, abs2(sx)); isnothing(ea) || return(something(ea))
    #!(Missing <: eltype(x)) && return(sem_cor(x, autocor_effective(x)))
    acfe = autocor_effective(x, ms; kwargs...)
    sem_cor(x, acfe, ms)
end,
function sem_cor(x, acfe, ms::MissingStrategy=PassMissing(); neff=nothing)
    sx = typeof(ms) <: HandleMissingStrategy ? std(skipmissing(x)) : std(x)::eltype(x)
    ea = early_var_return(x, abs2(sx)); isnothing(ea) || return(something(ea))
    #length(x) <= 1 && return(sx)
    if isnothing(neff); neff = effective_n_cor(x, acfe, ms); end
    σ2 = var_cor(x, acfe, ms; neff=neff)
    √(σ2/neff)
end

@doc raw"""
    var_cor(x, ms::MissingStrategy=PassMissing())
    var_cor(x, acf::AbstractVector, ms::MissingStrategy=PassMissing())

Estimate the variance for an autocorrelated series.

Zieba 2011 provide the following formula:
```math
Var(x) = \frac{n_{eff}}{n (n_{eff}-1)} \sum \left( x_i - \bar{x} \right)^2 
= {(n-1) n_{eff} \over n (n_{eff}-1)} Var_{uncor}(x)
```    

# Arguments
- `x`: An iterator of a series of observations
- `acf`: AutocorrelationFunction starting from lag 0. 
- `ms`: [`MissingStrategy`](@ref) passed to [`effective_n_cor`](@ref).
  Value of `SkipMissing()` speeds up computation compared to `ExactMissing()`,
  but leads to a negatively biased result with absolute value of the bias 
  increasing with the number of missings.

# Optional Arguments
- `neff`: may provide a precomputed number of observations for efficiency.
"""
function var_cor(x, ms::MissingStrategy=PassMissing(); neff=nothing) 
    varx = typeof(ms) <: HandleMissingStrategy ? var(skipmissing(x)) : var(x)::eltype(x)
    ea = early_var_return(x, varx); isnothing(ea) || return(something(ea))
    acf = autocor(x, ms)
    var_cor(x, autocor_effective(x, acf), ms)
end,
function var_cor(x, acfe, ms::MissingStrategy=PassMissing(); neff=nothing)
    varx = typeof(ms) <: HandleMissingStrategy ? var(skipmissing(x)) : var(x)::eltype(x)
    ea = early_var_return(x, varx); isnothing(ea) || return(something(ea))
    n = length(x)
    nmiss = count(ismissing.(x))
    nfin = n - nmiss
    if isnothing(neff); neff = effective_n_cor(x, acfe, ms); end
    σ2uncorr = var(skipmissing(x))
    # BLUE Var(x) for correlated: Zieba11 eq.(1) 
    σ2 = σ2uncorr*(nfin-1)*neff/(nfin*(neff-1))
end

function early_var_return(x, varx=var(x))
    ismissing(varx) && return(missing)
    !isfinite(varx) && return(varx)
    varx == zero(varx) && return(varx)
    nothing
end

@doc raw"""
    effective_n_cor(x, ms::MissingStrategy=ExactMissing()) 
    effective_n_cor(x, acf::AbstractVector, ms::MissingStrategy=ExactMissing())

Compute the number of effective observations for an autocorrelated series.

The formula in Zieba has been extended for missing values:
```math
n_{eff} = \frac{n_F}{1+{2 \over n_F} \sum_{k=1}^{min(n-1,n_k)} (n-k-m_k) \rho_k}
```
where ``n`` is the number of total records, ``n_F`` is the number of 
finite records, ``n_k`` is the nummber of components in the 
used autocorrelation function (``n-1`` if not estimated from the data)
,``\rho_k`` is the correlation, and 
``m_k`` is the number of pairs that contain a missing value for lag ``k``.

# Arguments
- `x`: An iterator of a series of observations
- `ms`: [`MissingStrategy`](@ref): set to `SkipMissing()` to speed up computation 
  (omitting [`count_forlags`](@ref) missing pairs) 
  but get a positively biased
  result with increasing bias with the number of missings. 
  This leads to a subsequent underestimated uncertainty of the sum or the mean.
- `acf`: AutocorrelationFunction starting from lag 0

# Examples
```jldoctest; output = false, setup = :(using LogNormals, Distributions, LinearAlgebra, Missings)
acf0 = [1,0.4,0.1]
Sigma = cormatrix_for_acf(100, acf0);
# 100 random variables each Normal(1,1)
dmn = MvNormal(ones(100), Symmetric(Sigma));
x = allowmissing(rand(dmn));    
x[11:20] .= missing
neff = effective_n_cor(x, acf0)
neff < 90
neff_biased = effective_n_cor(x, acf0, SkipMissing())
neff_biased > neff
# output
true
```
"""
function effective_n_cor(x, ms::MissingStrategy=ExactMissing())
    effective_n_cor(x, autocor(x,ms), ms)
end,
function effective_n_cor(x, acf::AbstractVector, ms::MissingStrategy=ExactMissing())
    # Zieba 2001 eq.(3)
    n = length(x)
    k = Base.OneTo(min(n,length(acf))-1) # acf starts with lag 0
    if ms == ExactMissing() && (Missing <: eltype(x))
        # see derivation in sem_cor.md
        # number of missing combinations due to missing in x
        #only julia 1.6 (m0, mk...) = count_forlags(ismissing,x,0:length(k))
        mka = count_forlags((x_i,x_iplusk)->ismissing(x_i) || ismissing(x_iplusk), x, 0:length(k))
        m0 = mka[1]
        mk = mka[2:end]
        nf = n - m0
        neff = nf/(1 + 2/nf*sum((n .- k .-mk) .* acf[k.+1]))  
    else
        neff = n/(1 + 2/n*sum((n .- k) .* acf[k.+1]))  
    end
end
