"""
    Helpers to fit LogNormal distributions.
"""
module LogNormals

export AbstractMoments, Moments, n_moments, moments,
    QuantilePoint, fit_mean_quantile, fit_mode_quantile,
    @qp, @qp_ll, @qp_l, @qp_m, @qp_u, @qp_uu, 
    @qs_cf90, @qs_cf95


using Distributions, StaticArrays

# general fitting to statistics
include("fitstats.jl")

# lognormal fitting
include("lognormal.jl")

# normal fitting
include("normal.jl")

# logitnormal fitting
include("logitnormal.jl")

# sum of lognormal random variables
include("sumlognormals.jl")

    
end # module
