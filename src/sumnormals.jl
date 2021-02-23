function Base.sum(dv::Union{Base.SkipMissing{DV},DV}; skipmissings::Val{B} = Val(false)) where 
    {DV <: AbstractDistributionVector{<:Normal}, B} 
    B == true && return(sum(skipmissing(dv)))
    # uncorrelated, only sum diagonal
    Ssum = s = zero(eltype(nonmissingtype(eltype(dv))))#zero(T)
    nterm = 0
    for d in dv
        μ,σ = params(d)
        Ssum += μ
        s += abs2(σ)
        nterm += 1
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    Normal(Ssum, √s)
end

# own method with argument isgapfilled, because now cannot use
# skipmissing any more and need to allocate nonmissing to 
function Base.sum(dv::AbstractDistributionVector{<:Normal}, 
    isgapfilled::AbstractVector{Bool}; 
    skipmissings::Val{B} = Val(false)) where B
    length(dv) == length(isgapfilled) || error(
        "argument gapfilled must have the same length as dv ($(length(dv))" *
        "but was $(length(isgapfilled)).")
    if B == true
        # need to allocate anyway with subsetting
        nonmissing = findall(.!ismissing.(dv))
        return(sum(
            @inbounds(dv[nonmissing]), @inbounds(isgapfilled[nonmissing])))
    end
    # uncorrelated, only sum diagonal
    Ssum = s = zero(eltype(nonmissingtype(eltype(dv))))
    nterm = 0
    for (i,d) in enumerate(dv)
        μ,σ = params(d)
        Ssum += μ
        if !isgapfilled[i]
            s += abs2(σ)
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    Normal(Ssum, √s)
end



