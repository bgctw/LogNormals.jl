"""
Tools that help using the LogNormal distribution.

- Fitting to various aggregate statistics
- Sum of correlated lognormal random variables
"""
module LogNormals

export AbstractMoments, Moments, n_moments, moments,
    QuantilePoint, 
    fit_mean_quantile, fit_mode_quantile, fit_median_quantile,
    @qp, @qp_ll, @qp_l, @qp_m, @qp_u, @qp_uu, 
    @qs_cf90, @qs_cf95,
    Σstar, σstar, 
    autocor_effective, effective_n_cor, count_for_lag, count_forlags, sem_cor, var_cor

import Random: rand, rand!
import Base: sum, size
import Statistics: mean
import StatsBase: coef, autocor

using StatsBase, Distributions, Missings, MissingStrategies
using BandedMatrices, MappedArrays, RecursiveArrayTools, FillArrays, OffsetArrays
using StaticArrays, LinearAlgebra
using Random
using SimpleTraits

# general fitting to statistics
include("fitstats.jl")

# lognormal fitting
include("lognormal.jl")

# normal fitting
include("normal.jl")

# logitnormal fitting
include("logitnormal.jl")

# standard error of the mean of correlated normal variables
include("semcor.jl")
    
end # module
