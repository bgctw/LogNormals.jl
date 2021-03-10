function sum_lognormals!(S, μ, σ, corr::AbstractMatrix, ms::MissingStrategy=PassMissing(); 
    corrlength = length(μ)-1)
    if typeof(ms) <: HandleMissingStrategy
        any(ismissing.(corr)) && error("cannot skip missings in correlation matrix")
        ismiss = ismissing.(μ) .| ismissing.(σ) 
        #.| vec(mapslices((x -> any(ismissing(x))), corr; dims = 1))
        if any(ismiss) 
            μ = @view μ[.!ismiss]
            σ = @view σ[.!ismiss]
            corr = @view corr[.!ismiss, .!ismiss]
        end           
    end
  ##details<< Implements estimation according to
  ## Messica A(2016) A simple low-computation-intensity model for approximating
  ## the distribution function of a sum of non-identical lognormals for
  ## financial applications. 10.1063/1.4964963
  #ρ
  nterm = length(μ) 
  @. S = exp(μ + abs2(σ)/2)
  Ssum = zero(μ[1])
  s = zero(μ[1])
  for i in 1:nterm
    Ssum += S[i]
    jstart = max(1, i - corrlength)
    jend = min(nterm, i + corrlength)
    for j in jstart:jend
        s += corr[i,j] * σ[i] * σ[j] * S[i] * S[j]
    end
  end
  σ2eff = s/abs2(Ssum)
  μ_sum = log(Ssum) - σ2eff/2
  LogNormal(μ_sum, √σ2eff)  
end

function sum_lognormals(μ, σ, corr::AbstractMatrix, ms::MissingStrategy=PassMissing(); 
    corrlength = length(μ)-1)
    sum_lognormals!(similar(μ), μ, σ, corr, ms, corrlength = corrlength)
end

function sum_lognormals(μ, σ, ms::MissingStrategy=PassMissing())
    if typeof(ms) <: HandleMissingStrategy
        ismiss = ismissing.(μ) .| ismissing.(σ)
        if any(ismiss) 
            μ = @view μ[.!ismiss]
            σ = @view σ[.!ismiss]
        end           
    end
    # uncorrelated, only sum diagonal
    nterm = length(μ)
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    Ssum = zero(μ[1])
    s = zero(μ[1])
    for i in 1:nterm
        Si = exp(μ[i] + abs2(σ[i])/2)
        Ssum += Si
        s += abs2(σ[i]) * abs2(Si)
    end
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end

function length_itr(x)
    typeof(Base.IteratorSize(x)) <: Union{Base.HasShape, Base.HasLength} && 
        return(length(x))
    count(x -> true, x)
end

# document and example skipmissing(ds)
function sum_lognormals(ds) 
    nterm = length_itr(ds)
    # uncorrelated, only sum diagonal
    nterm > 0 || error("Expected at least one nonmissing term, but mu = $μ")
    Ssum = s = zero(eltype(eltype(ds)))
    for i in 1:nterm
        μ,σ = params(ds[i])
        Si = exp(μ + abs2(σ)/2)
        Ssum += Si
        s += abs2(σ) * abs2(Si)
    end
    σ2eff = s/abs2(Ssum)
    μ_sum = log(Ssum) - σ2eff/2
    LogNormal(μ_sum, √σ2eff)
end
  
function sum_lognormals!(S, ds, corr::AbstractMatrix, ms::MissingStrategy=PassMissing(); 
    corrlength = length(ds)-1)
    nterm = length_itr(ds) 
    if typeof(ms) <: HandleMissingStrategy
        any(ismissing.(corr)) && error("cannot skip missings in correlation matrix")
        ifin = findall(.!ismissing.(ds))
        #.| vec(mapslices((x -> any(ismissing(x))), corr; dims = 1))
        if length(ifin) != nterm
            ds = @view ds[ifin]
            corr = @view corr[ifin, ifin]
        end           
    end
  ##details<< Implements estimation according to
  ## Messica A(2016) A simple low-computation-intensity model for approximating
  ## the distribution function of a sum of non-identical lognormals for
  ## financial applications. 10.1063/1.4964963
  #ρ
  nterm = length(μ) 
  @. S = exp(μ + abs2(σ)/2)
  Ssum = zero(μ[1])
  s = zero(μ[1])
  for i in 1:nterm
    Ssum += S[i]
    jstart = max(1, i - corrlength)
    jend = min(nterm, i + corrlength)
    for j in jstart:jend
        s += corr[i,j] * σ[i] * σ[j] * S[i] * S[j]
    end
  end
  σ2eff = s/abs2(Ssum)
  μ_sum = log(Ssum) - σ2eff/2
  LogNormal(μ_sum, √σ2eff)  
end

  
  