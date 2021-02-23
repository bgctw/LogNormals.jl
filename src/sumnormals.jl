function Base.sum(dv::Union{Base.SkipMissing{DV},DV}; 
    skipmissings::Val{B} = Val(false)) where 
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
    Ssum = s = Ssumnonfilled = zero(eltype(nonmissingtype(eltype(dv))))
    nterm = 0
    for (i,d) in enumerate(dv)
        μ,σ = params(d)
        Ssum += μ
        if !isgapfilled[i]
            Ssumnonfilled += μ
            s += abs2(σ)
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    relerr = √s/Ssumnonfilled
    Normal(Ssum, Ssum * relerr)
end

function Base.sum(dv::AbstractDistributionVector{D}, 
    acf::AbstractVector{<:Number}, 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal, DS<:eltype(D), SM}
    # currently use Matrix-method, Maybe implement fast with loop
    corrM = Symmetric(cormatrix_for_acf(length(dv), acf))
    return(sum_normals(
        dv, corrM, isgapfilled,skipmissings = skipmissings, storage = storage))
end


function Base.sum(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix}, 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal, DS<:eltype(D), SM}
    sum_normals(
        dv, corr, isgapfilled, storage = storage, skipmissings = skipmissings)
end

function sum_normals(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix}, 
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)); 
    storage::AbstractVector{Union{Missing,DS}} = 
        Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal,SM, DS<:eltype(D)}
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    # var_sum (s) is the sum across all Sigma, i.e. σT * corr * σ
    # missings and gapfilled values do not count -> set to zero
    if SM == true
        #@. storage = (x -> ismissing(x) ? zero(DS) : x)(σ)  
        @. storage = coalesce(σ, zero(DS))  
        Spure = disallowmissing(storage) 
    else
        Spure = disallowmissing(σ)
    end 
    Spure[isgapfilled] .= zero(DS)
    s = transpose(Spure) * corr * Spure
    Ssum = Ssumnonfilled = zero(DS)
    nterm = 0
    for (i,d) in enumerate(dv)
        ismissing(d) && continue
        μ,σ = params(d)
        Ssum += μ
        if !isgapfilled[i]
            Ssumnonfilled += μ
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    relerr = √s/Ssumnonfilled
    #@show Ssum, Ssumnonfilled, s, relerr
    Normal(Ssum, Ssum * relerr)
end



