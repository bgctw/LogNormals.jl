function count_forlag(pred,x,k::Integer)
    s = 0
    ax = OffsetArrays.no_offset_view(axes(x,1))
    for i in 1:(length(x)-k)
        if (pred(x[ax[i]])::Bool || pred(x[ax[i+k]])::Bool); s += 1; end
    end
    s
end

"""
    count_forlags(pred, x,lags)

Count the number of records in where both x and lagged version of x
fulfil a predicate.

# Arguments
- `x::Function(xi)::Bool`: The predicate to be applied to each record of x
- `x`: The series whose lags are inspected.
- `lags`: An iterator of Integer lag sizes
"""
count_forlags(pred, x,lags) = count_forlag.(tuple(pred), tuple(x),lags)


"""
    autocor(x::AbstractVector{Union{Missing,<:T}, lags; dmean, exactmissing=true}

Estimate autocorrelation function with presence of missing values.

# Arguments
- `exactmissing`: If set to `true` then devide the sum
    in the formula of the exepected value in the formula for the correlation
    at lag `k` by `n - nmissing` instead of `n`, 
    where `nimissing` is the number of records where there is a missing either
    in the original vector or its lagged version.

# Details    
    If `exactmissing=false` then the computation is faster, but
    it is more strongly biased low with increasing number of missings. 
    Note that `StatsBase.autocor` uses devision by `n` instead of 'n-k', the
    true length of the vectors correlated at lag `k` resulting in 
    low-biased correlations of higher lags for numerical stability reasons.
"""
function autocor(x::AbstractVector{Union{Missing,T}}, 
    lags::AbstractVector{<:Integer}=StatsBase.default_autolags(size(x,1)); 
    demean::Bool=true, exactmissing::Bool=true, kwargs...) where T <: Real
    # only set to zero after demeaning, otherwise correlation estimates increase
    z::Vector{Union{Missing,T}} = demean ? x .- mean(skipmissing(x)) : x
    # replace missing by zero: new type will not match signature of current function 
    zpure = coalesce.(z, zero(z))::Vector{T} 
    acf = autocor(zpure,lags;demean=false,kwargs...)
    !exactmissing && return(acf)
    # correct for sum has been devided by a larger number of terms
    # (including terms of missing/0) 
    lx = length(x)
    #zcorr = nterm_forlag(x,0)/lx
    zcorr_can = lx - count_forlag(ismissing,x,0)#nterm_forlag(x,0)
    for (i,k) in enumerate(lags)
        #autocor code devides by lx instead of (lx-k)
        #https://github.com/JuliaStats/StatsBase.jl/issues/273#issuecomment-307560660
        # this cancels with lx of zcorr
        #acf[i] *= (lx - k)/nterm_forlag(x,k) * zcorr
        acf[i] *= zcorr_can/(lx - count_forlag(ismissing,x,k)) 
    end
    acf
end

function autocor_effective(x, acf = autocor(x))
    #maybe implement a more efficient version that computes only the
    # first lags and further lags if not found negative correlation
    i = findfirst(x -> x <= 0.0, acf)
    isnothing(i) && return(acf)
    acf[1:(i-1)]
end

function sem_cor(x; exactmissing::Bool=true) 
    Missing <: eltype(x) || return(sem_cor(x, autocor_effective(x)))
    acfe = autocor_effective(x, autocor(x, exactmissing=exactmissing))
    sem_cor(x, acfe; exactmissing=exactmissing)
end

function sem_cor(x, acfe; exactmissing::Bool=true)
    n = length(x)
    n <= 1 && return(std(x))
    nmiss = count(ismissing.(x))
    nfin = n - nmiss
    neff = effective_n_cor(x, acfe; exactmissing=exactmissing)
    σ2uncorr = var(skipmissing(x))
    # BLUE Var(x) for correlated: Zieba11 eq.(1) 
    #σ2 = σ2uncorr*(nfin-1)*neff/(nfin*(neff-1))
    #√(σ2/neff)
    √(σ2uncorr*(nfin-1)/(nfin*(neff-1))) # one neff cancelled
end

function effective_n_cor(x, acf::AbstractVector; exactmissing::Bool=true)
    # Zieba 2001 eq.(3)
    n = length(x)
    k = Base.OneTo(min(n,length(acf))-1) # acf starts with lag 0
    if exactmissing && (Missing <: eltype(x))
        # see derivation in sem_cor.md
        # number of missing combinations due to missing in x
        #only julia 1.6 (m0, mk...) = count_forlags(ismissing,x,0:length(k))
        mka = count_forlags(ismissing,x,0:length(k))
        m0 = mka[1]
        mk = mka[2:end]
        nf = n - m0
        neff = nf/(1 + 2/nf*sum((n .- k .-mk) .* acf[k.+1]))  
    else
        neff = n/(1 + 2/n*sum((n .- k) .* acf[k.+1]))  
    end
end
