# Sum of correlated lognormal random variables


<!-- ```@meta
function boot_dvsums_acf(dv, acf, nboot = 10_000)
    # test sum formula by bootstrap sample
    μ, σ = params(dv)
    Sigma = Diagonal(σ) * cormatrix_for_acf(length(dv), acf) * Diagonal(σ);
    dn = MvNormal(disallowmissing(μ), Symmetric(Sigma));
    x = rand(dn, nboot) .|> exp
    sums = vec(sum(x, dims = 1))
    #density(sums)
    drsum = fit(LogNormal, sums)
end
```

```jldoctest; output = false, setup = :(using Distributions,LogNormals,StatsPlots,Plots)
    mu = log.([110,100,80,120,160.0])
    sigma = log.([1.2,1.5,1.1,1.3,1.1])
    acf1 = [0.4,0.1]
    dv = SimpleDistributionVector(LogNormal{eltype(mu)}, mu, sigma);
    drsum = boot_dvsums_acf(dv, acf1) # boot sum across random numbers

    dsum = sum(dv, acf1)

    @test dsum ≈ drsum rtol = 0.02
``` -->

```@docs
sum(Union{Base.SkipMissing{DV},DV}; skipmissings::Val{B} = Val(false)) where 
    {T, DV <: AbstractDistributionVector{LogNormal{T}}, B}
```

