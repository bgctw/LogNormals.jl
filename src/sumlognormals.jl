function sum(dv::AbstractDistributionVector{<:LogNormal}; 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)),
    skipmissings::Val{B} = Val(false)) where B
    length(dv) == length(isgapfilled) || error(
        "argument gapfilled must have the same length as dv ($(length(dv))" *
        "but was $(length(isgapfilled)).")
    if B == true
        nonmissing = findall(.!ismissing.(dv))
        !isempty(nonmissing) && return(sum(
            @inbounds(dv[nonmissing]), 
            isgapfilled=@inbounds(isgapfilled[nonmissing])))
    end
    # uncorrelated, only sum diagonal
    Ssum = s = Ssumnonfilled = zero(eltype(nonmissingtype(eltype(dv))))#zero(T)
    nterm = 0
    for (i,d) in enumerate(dv)
        μ,σ = params(d)
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        if !isgapfilled[i]
            Ssumnonfilled += Si
            s += abs2(σ) * abs2(Si)
            nterm += 1
        end
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    σ2eff = s/abs2(Ssumnonfilled)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end


function sum(dv::AbstractDistributionVector{D}, 
    acf::AutoCorrelationFunction; 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)), 
    storage::AbstractVector{Union{Missing,ST}} = 
        Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false),    
    method::Val{M} = Val(:vector)) where 
    {D<:LogNormal, SM, ST<:eltype(D), M} 
    #storage = Vector{Union{Missing,eltype(D)}}(undef, length(dv))
    if M == :vector
        return(sum_lognormals(
            dv, acf, isgapfilled=isgapfilled, storage = storage, 
            skipmissings = skipmissings))
    end
    if M == :bandedmatrix
        corrM = Symmetric(cormatrix_for_acf(length(dv), coef(acf)))
        return(sum_lognormals(
            dv, corrM, isgapfilled=isgapfilled,skipmissings = skipmissings, 
            storage = storage))
    end
    error("Unknown method $method")
end

function sum_lognormals(dv::AbstractDistributionVector{D}, 
    acf::AutoCorrelationFunction; 
    isgapfilled::AbstractVector{Bool} = Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SK} = Val(false)) where 
    #{D<:LogNormal, ST<:, B}
    {D<:LogNormal, DS<:eltype(D), SK}
    #details<< Implements estimation according to
    # Messica A(2016) A simple low-computation-intensity model for approximating
    # the distribution function of a sum of non-identical lognormals for
    # financial applications. 10.1063/1.4964963
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    coef_acf = coef(acf)
    corrlength = length(coef_acf)
    acfm = vcat(reverse(coef_acf), 1, coef_acf)
    n = length(μ)
    @. storage = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing.(storage))
    skipmissings != Val(true) && nmissing != 0 && error(
        "Found missing values. Use argument 'skipmissings = Val(true)' " *
        "to sum over nonmissing.")
    # 0 in storage has the effect of not contributing to Ssum nor s 
    # For excluding terms of gapfilled records, must make 
    # to either set storage[gapfilled] to zero or check in the product
    Spure = disallowmissing(replace(storage, missing => 0.0))
    σpure = disallowmissing(replace(σ, missing => 0.0))
    Ssum = sum(Spure)
    # after computing Ssum, can also set gapfilled to zero
    Spure[isgapfilled] .= 0.0
    s = Ssumunfilled = zero(eltype(σpure))
    for i in 1:n
        iszero(Spure[i]) && continue # nothing added
        Ssumunfilled += Spure[i]
        jstart = max(1, i - corrlength)
        jend = min(n, i + corrlength)
        for j in jstart:jend
            acf_ind = (j-i + corrlength +1)
            # sij will be zero if sigma or storage is missing (replaced by zero)
            # Sj moved to start to help multiplication by early zero
            sij = Spure[j] * acfm[acf_ind] * σpure[i] * σpure[j] * Spure[i] 
            s += sij
        end
    end
    σ2eff = s/abs2(Ssumunfilled)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, n - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

function sum(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix}; 
    isgapfilled::AbstractArray{Bool,1}=Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
       Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:LogNormal, DS<:eltype(D), SM}
    sum_lognormals(
        dv, corr, isgapfilled=isgapfilled, storage = storage, 
        skipmissings = skipmissings)
end

function sum_lognormals(dv::AbstractDistributionVector{D}, 
    corr::Symmetric{DS,<:AbstractMatrix};
    isgapfilled::AbstractArray{Bool,1} = Falses(length(dv)),
    storage::AbstractVector{Union{Missing,DS}} = 
        Vector{Union{Missing,eltype(D)}}(undef, length(dv)),
    skipmissings::Val{SM} = Val(false)) where 
    {D<:LogNormal,SM, DS<:eltype(D)}
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    # storage = allowmissing(similar(μ))
    @. storage = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing, storage)
    anymissing = nmissing != 0
    SM != true && anymissing && error(
         "Found missing values. Use argument 'skipmissings = Val(true)' " *
         "to sum over nonmissing.")
    Ssum::nonmissingtype(eltype(storage)) = sum(skipmissing(storage))
    # gapfilled records only used for Ssum, can set the to 0 now
    # so they do not contribute to s and Ssumfilled for computation of σ2eff
    storage[isgapfilled] .= 0
    Ssumunfilled::nonmissingtype(eltype(storage)) = sum(skipmissing(storage))
    @. storage = σ * storage  # do only after Ssum
    # setting storage to zero results in summing zero for missing records
    # which is the same as filtering both storage and corr
    anymissing && replace!(storage, missing => 0.0)
    #s = transpose(disallowmissing(storage)) * corr * disallowmissing(storage)
    #Spure = view_nonmissing(storage) # non-allocating
    Spure = disallowmissing(storage) # allocating - tested: is faster than the view
    s = transpose(Spure) * corr * Spure
    σ2eff = s/abs2(Ssumunfilled)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, length(storage) - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

