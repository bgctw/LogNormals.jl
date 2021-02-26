function sum(dv::AbstractDistributionVector{<:Normal}; 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)), 
    skipmissings::Val{B} = Val(false)) where B
    length(dv) == length(isgapfilled) || error(
        "argument gapfilled must have the same length as dv ($(length(dv))" *
        "but was $(length(isgapfilled)).")
    if B == true
        # need to allocate anyway with subsetting
        nonmissing = findall(.!ismissing.(dv))
        if !isempty(nonmissing) 
            return(sum(@inbounds(dv[nonmissing]), 
                isgapfilled = @inbounds(isgapfilled[nonmissing])))
        end
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

function sum(dv::AbstractDistributionVector{D}, 
    acf::AutoCorrelationFunction;
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal, DS<:eltype(D), SM}
    # currently use Matrix-method, Maybe implement fast with loop
    corrM = Symmetric(cormatrix_for_acf(length(dv), coef(acf)))
    return(sum_normals(
        dv, corrM, 
        isgapfilled=isgapfilled,skipmissings = skipmissings, storage = storage))
end


function sum(dv::AbstractDistributionVector{D}, 
    corr::Symmetric; 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:Normal, DS<:eltype(D), SM}
    sum_normals(
        dv, corr, 
        isgapfilled=isgapfilled, storage=storage, skipmissings=skipmissings)
end

function sum_normals(dv::AbstractDistributionVector{D}, 
    corr::Symmetric; 
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)),
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
    Normal(Ssum, Ssum * relerr)
end


mean(dv::AbstractDistributionVector{<:Normal}; kwargs...) =
    mean_normals(dv; kwargs...)
mean(dv::AbstractDistributionVector{<:Normal}, corr::Symmetric; kwargs...) =
    mean_normals(dv, corr; kwargs...)
mean(dv::AbstractDistributionVector{<:Normal}, acf::AutoCorrelationFunction; 
    kwargs...) =
    mean_normals(dv, acf; kwargs...)

function mean_normals(dv::AbstractDistributionVector{<:Normal}, x...; 
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)), kwargs...) 
    ds = sum(dv, isgapfilled=isgapfilled, x...; kwargs...)
    n = count(x -> !ismissing(x), dv)
    nobs = count((t->!ismissing(t[1]) && !t[2]), zip(dv, isgapfilled))
    # relative error scales by 1/√nobs
    Normal(ds.μ/n, ds.σ/(n*√nobs))
end


