

# function length_itr(x)
#     typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
#         return(length(x))
#     count(x -> true, x)
# end

function Base.sum(dv::DSM; skipmissings::Val{B} = Val(false)) where 
    DSM <: Union{Base.SkipMissing{DV},DV} where 
    {DV <: AbstractDistributionVector{LogNormal{T}}, B} where T
    B == true && return(sum(skipmissing(dv)))
    # uncorrelated, only sum diagonal
    Ssum = s = zero(T)
    nterm = 0
    for d in dv
        μ,σ = params(d)
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        s += abs2(σ) * abs2(Si)
        nterm += 1
    end
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end

function cormatrix_for_acf(n::Int,acf::AbstractVector) 
    nacf::Int = length(acf)
    corrM = BandedMatrix{Float64}(undef, (n,n), (nacf,nacf))
    corrM[band(0)] .= 1
    for i in 1:nacf
      corrM[band(i)] .= corrM[band(-i)] .= acf[i]
    end
    corrM
end

function Base.sum(dv::AbstractDistributionVector{<:LogNormal}, 
    acf::AbstractVector; 
    skipmissings::Val{B} = Val(false), method::Val{S} = Val(:vector)) where 
    {B, S} 
    storage = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(dv))
    if method == Val(:vector) 
        return(sum_lognormals!(storage, dv, acf, skipmissings = skipmissings))
    end
    if method == Val(:bandedmatrix)
        corrM = cormatrix_for_acf(length(dv), acf)
        return(sum_lognormals!(storage, dv, corrM, skipmissings = skipmissings))
    end
    error("Unknown method $method")
end

function sum_lognormals!(S::Vector{Union{Missing,T}}, dv::DV, 
    acf::AbstractVector; 
    skipmissings::Val{B} = Val(false)) where 
    {DV <: AbstractDistributionVector{D}, B} where D<:LogNormal where T
    #details<< Implements estimation according to
    # Messica A(2016) A simple low-computation-intensity model for approximating
    # the distribution function of a sum of non-identical lognormals for
    # financial applications. 10.1063/1.4964963
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    corrlength = length(acf)
    acfm = vcat(reverse(acf), 1, acf)
    n = length(μ)
    @. S = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing.(S))
    skipmissings != Val(true) && nmissing != 0 && error(
        "Found missing values. Use argument 'skipmissings = Val(true)' " *
        "to sum over nonmissing.")
    Spure = disallowmissing(replace(S, missing => 0.0))
    σpure = disallowmissing(replace(σ, missing => 0.0))
    Ssum = s = zero(eltype(σpure))
    for i in 1:n
        Spure[i] == zero(Spure) && continue # nothing added
        Ssum += Spure[i]
        jstart = max(1, i - corrlength)
        jend = min(n, i + corrlength)
        for j in jstart:jend
            acf_ind = (j-i + corrlength +1)
            # sij will be zero if sigma or S is missing (replaced by zero)
            # Spure moved to start because may be zero 
            sij = Spure[j] * acfm[acf_ind] * σpure[i] * σpure[j] * Spure[i] 
            s += sij
        end
    end
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, n - nmissing
    LogNormal(μ_sum, √σ2eff)  
end


function Base.sum(dv::AbstractDistributionVector{<:LogNormal}, 
    corr::AbstractMatrix; skipmissings::Val{B} = Val(false)) where B
    S = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(dv))
    sum_lognormals!(S, dv, corr, skipmissings = skipmissings)
end

#view_nonmissing(x) = of_eltype(nonmissingtype(eltype(x)),x)

function sum_lognormals!(S, dv, corr::AbstractMatrix; 
    skipmissings::Val{l} = Val(false)) where l
    μ = params(dv, Val(1))
    σ = params(dv, Val(2))
    # S = allowmissing(similar(μ))
    @. S = exp(μ + abs2(σ)/2)
    nmissing = count(ismissing, S)
    anymissing = nmissing != 0
    l != true && anymissing && error(
         "Found missing values. Use argument 'skipmissings = Val(true)' " *
         "to sum over nonmissing.")
    Ssum::nonmissingtype(eltype(S)) = sum(skipmissing(S))
    @. S = σ * S  # do only after Ssum
    # setting S to zero results in summing zero for missing records
    # which is the same as filtering both S and corr
    anymissing && replace!(S, missing => 0.0)
    #s = transpose(disallowmissing(S)) * corr * disallowmissing(S)
    #Spure = view_nonmissing(S) # non-allocating
    Spure = disallowmissing(S) # allocating - tested: is faster than the view
    s = transpose(Spure) * corr * Spure
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    #@show Ssum, s, length(S) - nmissing
    LogNormal(μ_sum, √σ2eff)  
end

