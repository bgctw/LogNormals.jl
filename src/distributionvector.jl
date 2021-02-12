## AbstractDistributionVector
"""
    AbstractDistributionVector{D <: Distribution}

Is any type able iterating the same type of Distribution.
Parametrized by a `Distribution` defining the type of the Distribution.

Items may be missing. Hence the element type of the Iterator is `Union{Missing,D}`.

Specific implementations, such as [DistributionVector](@ref) need
to implement methods `length` and `getindex`.

Method `params(dv, i)` returns the AbstractVector of the ith parameter of each
distribution. The default implementation for the abstract type, 
is based on the `params` method applied to each the distribution. 
Concrete Types should should provide a more efficient 
specialization.
"""
abstract type AbstractDistributionVector{D <: Distribution} end

Base.eltype(::Type{<:AbstractDistributionVector{D}}) where D = Union{Missing,D}

function Base.iterate(dv::AbstractDistributionVector, state=1) 
    state > length(dv) ? nothing : (dv[state], state+1)
end

function Base.iterate(rds::Iterators.Reverse{AbstractDistributionVector}, 
    state=length(rds.itr))  
    state < 1 ? nothing : (rds.itr[state], state-1)
end

StatsBase.params(dv::AbstractDistributionVector, i::Integer) = getindex.(params.(dv),i)



## DistributionVector   
"""
    DistributionVector{D <: Distribution, V}

Is an Vector-of-Tuples based implementation of [AbstractDistributionVector](@ref).
It is parametrized by the type of `Distribution` D, the type of abstract vector
of tuples V.

Method `params` returns a (Npar x N) AbstractArray of distribution parameters
"""
struct DistributionVector{D <: Distribution, V <: AbstractVector} <: 
    AbstractDistributionVector{D} 
    params::V
    # inner constructor checking ?
end

function DistributionVector(::Type{D}, params::V) where 
{D<:Distribution,  V<:AbstractVector} 
    Missing <: eltype(V) || error(
        "Expected type of parameters to allow for missing." *
        " Can you use 'allowmissing' in constructing the DistributionVector?")
    eltype(V) <: Union{Missing, <:Tuple} || error(
        "Expected type of parameters of 'Union{Missing, <:Tuple}' "*
        " but got $(eltype(V)).")
    DistributionVector{D, V}(params)
end

function DistributionVector(dv::Vararg{Union{Missing,D},N}) where {D<:Distribution, N} 
    N == 0 && error(
        "Provide at least one distribution in DistributionVector(x...).")
    parms = collect(passmissing(params).(dv))
    DistributionVector(D, allowmissing(parms))
end

function DistributionVector(::Type{D}, pvec::Vararg{Any,N}) where {D<:Distribution, N} 
    # if one parameter has missing, the entire tuple must be set to missing
    replacemissing(x) = (ismissing(x) || any(ismissing.(x))) ? missing : x
    parms = collect(replacemissing.(zip(pvec...)))
    DistributionVector(D, allowmissing(parms))        
end


Base.length(dv::DistributionVector) = length(dv.params)

function Base.getindex(dv::DistributionVector{D,V},i::Int)::Union{Missing, D} where {D,V}
    params_i = dv.params[i]
    ismissing(params_i) && return missing
    D(params_i...)
end

function StatsBase.params(dv::DistributionVector, i::Integer) 
   # mappedarray(e -> passmissing(getindex)(e,i), dv.params)
   # currentl does not work, see https://github.com/JuliaArrays/MappedArrays.jl/issues/40
   passmissing(getindex).(dv.params,i)
end

