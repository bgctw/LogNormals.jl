

# function length_itr(x)
#     typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
#         return(length(x))
#     count(x -> true, x)
# end

nparams(::Type{<:LogNormal}) = 2
nparams(::Type{<:Normal}) = 2
nparams(::Type{<:LogitNormal}) = 2



# function Base.sum(ds::DSM; skipmissings::Val{B} = Val(false)) where 
#     DSM <: Union{Base.SkipMissing{DS},DS} where 
#     {DS <: AbstractDistributionSequence{LogNormal{T}}, B} where T
#     skipmissings == Val(true) && return(sum(skipmissing(ds)))
#     # uncorrelated, only sum diagonal
#     Ssum = s = zero(T)
#     nterm = 0
#     for d in ds
#         μ,σ = params(d)
#         Si = exp(μ + abs2(σ)/2)
#         Ssum += Si
#         s += abs2(σ) * abs2(Si)
#         nterm += 1
#     end
#     nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
#     σ2eff = s/abs2(Ssum)
#     μ_sum = log(Ssum) - σ2eff/2
#     LogNormal(μ_sum, √σ2eff)
# end

# function cormatrix_for_acf(n::Int,acf::AbstractVector) 
#     nacf::Int = length(acf)
#     corrM = BandedMatrix{Float64}(undef, (n,n), (nacf,nacf))
#     corrM[band(0)] .= 1
#     for i in 1:nacf
#       corrM[band(i)] .= corrM[band(-i)] .= acf[i]
#     end
#     corrM
# end

# function Base.sum(ds::DS, acf::AbstractVector; skipmissings::Val{B} = Val(false), method::Val{S} = Val(:vector)) where 
#     {DS <: AbstractDistributionSequence{LogNormal}, B, S} 
#     storage = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(ds))
#     if method == Val(:vector) 
#         return(sum_lognormals!(storage, ds, acf, skipmissings = skipmissings))
#     end
#     if method == Val(:bandedmatrix)
#         corrM = cormatrix_for_acf(length(ds), acf)
#         return(sum_lognormals!(storage, ds, corrM, skipmissings = skipmissings))
#     end
#     error("Unknown method $method")
# end


# function sum_lognormals!(S::Vector{Union{Missing,T}}, ds::DS, acf::AbstractVector; 
#     skipmissings::Val{B} = Val(false)) where 
#     {DS <: AbstractDistributionSequence{D}, B} where D<:LogNormal where T
#     ##details<< Implements estimation according to
#     ## Messica A(2016) A simple low-computation-intensity model for approximating
#     ## the distribution function of a sum of non-identical lognormals for
#     ## financial applications. 10.1063/1.4964963
#     parms = params(ds)
#     μ = @view parms[1,:]
#     σ = @view parms[2,:]
#     corrlength = length(acf)
#     acfm = vcat(reverse(acf), 1, acf)
#     n = size(parms,2)
#     @. S = exp(μ + abs2(σ)/2)
#     nmissing = count(ismissing.(S))
#     S2 = disallowmissing(S)
#     S2[1] = 3.0
#     replace!(S, missing => 0.0)
#     Ssum = s = 0.0 #zero(eltype(LogNormal))
#     for i in 1:n
#         Ssum += S[i]
#         jstart = max(1, i - corrlength)
#         jend = min(n, i + corrlength)
#         for j in jstart:jend
#             acf_ind = (j-i + corrlength +1)
#             sij = acfm[acf_ind] * σ[i] * σ[j] * S[i] * S[j]
#             if !ismissing(sij) 
#                 s += sij
#             end
#         end
#     end
#     skipmissings != Val(true) && nmissing != 0 && error(
#         "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
#     #σ2eff::eltype(LogNormal) = s/abs2(Ssum)
#     #μ_sum::eltype(LogNormal) = log(Ssum) - σ2eff/2
#     σ2eff = s/abs2(Ssum)
#     μ_sum = log(Ssum) - σ2eff/2
#     #@show Ssum, s, n - nmissing
#     LogNormal(μ_sum, √σ2eff)  
# end


# function Base.sum(ds::DS, corr::AbstractMatrix; skipmissings::Val{B} = Val(false)) where 
#     {DS <: AbstractDistributionSequence{LogNormal}, B} 
#     S = Vector{Union{Missing,eltype(LogNormal)}}(undef, length(ds))
#     sum_lognormals!(S, ds, corr, skipmissings = skipmissings)
# end

# view_nonmissing(x) = of_eltype(nonmissingtype(eltype(x)),x)

# function sum_lognormals!(S, ds, corr::AbstractMatrix; 
#     skipmissings::Val{l} = Val(false)) where l
#     parms = params(ds)
#     μ = @view parms[1,:]
#     σ = @view parms[2,:]
#     # S = allowmissing(similar(μ))
#     @. S = exp(μ + abs2(σ)/2)
#     nmissing = count(ismissing, S)
#     anymissing = nmissing != 0
#     skipmissings != Val(true) && anymissing && error(
#          "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
#     Ssum::nonmissingtype(eltype(S)) = sum(skipmissing(S))
#     @. S = σ * S  # do only after Ssum
#     # setting S to zero results in summing zero for missing records
#     # which is the same as filtering both S and corr
#     anymissing && replace!(S, missing => 0.0)
#     #s = transpose(disallowmissing(S)) * corr * disallowmissing(S)
#     #Spure = view_nonmissing(S) # non-allocating
#     Spure = disallowmissing(S) # allocating - tested: is faster than the view
#     s = transpose(Spure) * corr * Spure
#     σ2eff = s/abs2(Ssum)
#     μ_sum = log(Ssum) - σ2eff/2
#     #@show Ssum, s, length(S) - nmissing
#     LogNormal(μ_sum, √σ2eff)  
# end

