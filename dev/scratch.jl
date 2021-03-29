using LogNormals
using SimpleTraits
using MissingStrategies

macro handlemissings_stub2(
    # omit the forwarding method
    fun,pos_missing=1, pos_strategy=pos_missing +1, type_missing=Any,
    defaultstrategy=nothing, argname_strategy=:ms, suffix = "_hm",
    )
    #gens = (mgen.forwarder, mgen.missingstrategy_notsuperofeltype)
    gens = (mgen.missingstrategy_notsuperofeltype,)
    #gens = (mgen.forwarder, )
    #gens = ()
    esc(handlemissings(fun, pos_missing, pos_strategy, type_missing,
        defaultstrategy, gens, argname_strategy, suffix))
 end
# #tmp2 = @macroexpand 
# macro handlemissings_stub2(
#     fun,pos_missing=1, pos_strategy=pos_missing +1, type_missing=Any,
#     defaultstrategy=nothing,suffix = "_hm",
#     )
#     #gens = (mgen.forwarder, mgen.missingstrategy_notsuperofeltype)
#     gens = (mgen.missingstrategy_notsuperofeltype,)
#     esc(handlemissings(fun, pos_missing, pos_strategy, type_missing,
#         defaultstrategy, gens, suffix))
#   end
# @handlemissings_stub2(
#     autocor(x::AbstractVector{<:Real}; demean::Bool=true) = 0,
#     1,2,AbstractVector{<:Union{Missing,Real}}, nothing, 
# )
@handlemissings_stub2(
    autocor(x::AbstractVector{<:Real}, lags::AbstractVector{<:Integer}; demean::Bool=true) = 0,
    1,2,AbstractVector{<:Union{Missing,Real}}, nothing, 
)
#autocor(x::AbstractVector{<:Union{Missing, Real}}, ms::MissingStrategies.MissingStrategy; kwargs...) = 1


# @handlemissings_stub(
#     autocor(x::AbstractVector{<:Real}; demean::Bool=true) = 0,
#     1,3,AbstractVector{<:Union{Missing,Real}}, PassMissing(),
# )

  

# tmp3 = @macroexpand @traitfn function autocor_hmnolag(ms::MissingStrategies.MissingStrategy, x::::IsEltypeSuperOfMissing; 
#     demean::Bool=true
#     )
#     "Super of missing defined in LogNormals"
#     # lags not given -> forward to method with standard lag size
#     #autocor_hmlag(ms, x, StatsBase.default_autolags(size(x,1)); demean)
# end
# @traitfn function autocor_hm(ms::MissingStrategy, x::::IsEltypeSuperOfMissing; 
#     demean::Bool=true
#     )
#     "Super of missing defined in LogNormals"
#     # lags not given -> forward to method with standard lag size
#     #autocor_hmlag(ms, x, StatsBase.default_autolags(size(x,1)); demean)
# end
@traitfn function autocor_hm(ms::MissingStrategy, x::::IsEltypeSuperOfMissing, 
    lags::AbstractVector{<:Integer}; 
    demean::Bool=true
    )
    "Super of missing defined in LogNormals"
    # lags not given -> forward to method with standard lag size
    #autocor_hmlag(ms, x, StatsBase.default_autolags(size(x,1)); demean)
end

# module Test1
#     using SimpleTraits, MissingStrategies
#     @traitfn function autocor_hmnolag(ms::MissingStrategies.MissingStrategy, x::::IsEltypeSuperOfMissing; 
#         demean::Bool=true
#         )
#         "Super of missing defined in Test1"
#         # lags not given -> forward to method with standard lag size
#         #autocor_hmlag(ms, x, StatsBase.default_autolags(size(x,1)); demean)
#     end
# end

# MissingStrategies.@testdefinetrait_esc_opt()
# @traitfn function f1_esc_opt(x::::!(SimpleTraits.BaseTraits.IsBits); demean::Bool=true)
#     "f1_esc_opt !IsBist from LogNormals"  
# end

using Missings
m = allowmissing(rand(Float64,(3,3)))
convert.(Float64,m) # convert preserves dimensions

