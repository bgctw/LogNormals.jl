using LogNormals
using Test, Distributions, LinearAlgebra, StatsBase, Missings, Random
using OffsetArrays, RecursiveArrayTools

function boot_sem_cor()
    #using FillArrays
    ##using BandedMatrices
    ##using Bootstrap
    nsum = 100
    probgap = 0.6
    acf0 = [1,0.4,0.1]
    Sigma = cormatrix_for_acf(nsum, acf0);
    dmn = MvNormal(ones(nsum), Symmetric(Sigma));
    #always use the same gaps
    igap = sample(1:nsum,trunc(Int,nsum*probgap), replace=false)
    b = allowmissing(rand(dmn));
    b[igap] .= missing;
    sem_cor(b, acf0, ExactMissing())
    nboot = 10_000
    resboot = map(1:nboot) do iboot
        b = allowmissing(rand(dmn));
        b[igap] .= missing
        sum(skipmissing(b)), mean(skipmissing(b)), sem_cor(b, acf0, ExactMissing())
    end;
    sums_boot, means_boot, semcor_boot = vectuptotupvec(resboot);
    mean(sums_boot), std(sums_boot)
    mean(means_boot), std(means_boot)
    dvn = SimpleDistributionVector(Fill(Normal(1,1), nsum)...);
    dvn[ismissing.(b)] .= missing;
    dsum = sum(dvn, AutoCorrelationFunction(acf0), SkipMissing())
    dmean = mean(dvn, AutoCorrelationFunction(acf0), SkipMissing())
    @test dmean.σ ≈ std(means_boot) atol=0.01
    # mean over Normals with missings works
    function sem_cor_benchmark(x, acfe, ms::MissingStrategy=ExactMissing())
        n = length(x)
        n <= 1 && return(std(x))
        nmiss = count(ismissing.(x))
        nfin = n - nmiss
        neff = effective_n_cor(x, acfe)
        σ2uncorr = var(skipmissing(x))
        # BLUE Var(x) for correlated: Zieba11 eq.(1) 
        σ2 = σ2uncorr*(nfin-1)*neff/(nfin*(neff-1))
        if ms === ExactMissing() && (nmiss != 0)
            σ = Fill(√σ2, n)
            dv = ParamDistributionVector(Normal{nonmissingtype(eltype(x))}, x, σ)
            acfes = n < length(acfe) ? view(acfe,1:n) : acfe
            dm = mean(dv, AutoCorrelationFunction(acfes), SkipMissing())
            #@show acfes, dm, length(dv), count(ismissing,dv), count(ismissing,x)
            return(dm.σ)
        else
            return(√(σ2/neff))
        end
    end
    res_sem_cor = sem_cor(b, acf0)
    @test res_sem_cor ≈ std(means_boot) atol=0.02
    res_sem_cor2a = sem_cor_benchmark(b, acf0, ms=ExactMissing())
    @test res_sem_cor2a ≈ std(means_boot) atol=0.02
    res_sem_cor2b = sem_cor_benchmark(b, acf0, ms=SkipMissing())
    @test res_sem_cor2b ≈ std(means_boot) atol=0.02
    #using BenchmarkTools: @btime
    @btime sem_cor_benchmark($b, $acf0, ms=ExactMissing())
    @btime sem_cor_benchmark($b, $acf0, ms=SkipMissing())
    # neff based method much faster
    function tmpf()
        #using StatsPlots
        plot(density(means_boot))
        vline!([mean(means_boot), quantile(means_boot,[0.025,0.975])...])
        #
        # plot(density(means_boot), label = "density mean")
        # vline!([mean(means_boot), quantile(means_boot,[0.025,0.975])...], label = "stats mean")
        density(semcor_boot, label = "density sem_cor")
        vline!([mean(semcor_boot), quantile(semcor_boot,[0.025,0.975])...], label = "stats sem_cor")
        vline!([std(means_boot)], label = "std(means_bootstrap)") # true mu and sigma
        vline!([dmean.σ], label = "sum_normals - true distribs") # true mu and sigma
    end
end

function inspect_lags_autocorrelation()
    nsum = 1_000
    acf0 = [1,0.4,0.1]
    Sigma = cormatrix_for_acf(nsum, acf0);
    dmn = MvNormal(ones(nsum), Symmetric(Sigma));
    a = rand(dmn);
    probgap = 0.4
    clustersize = 2
    #clustersize = 50
    nboot = 1_000
    resboot = map(1:nboot) do iboot
        igap0 = sample(0:(nsum-clustersize),trunc(Int,nsum*probgap/clustersize), replace=true);
        b = allowmissing(a);
        # outer product see https://stackoverflow.com/a/44592419
        #igap0.+Base.OneTo(clustersize)' 
        b[igap0.+Base.OneTo(clustersize)'] .= missing;
        ae = autocor(b, 1:2, ms=ExactMissing())
        af = autocor(b, 1:2, ms=SkipMissing())
        (ae, af)
    end;
    ae, af = VectorOfArray.(vectuptotupvec(resboot));
    mean(ae; dims=2)
    mean(af; dims=2)
    autocor(a, 1:2)

end

function benchmark_count_forlags()
    #using BenchmarkTools: @btime
    function count_forlags_slower(pred, x,lags)
        cnt = zeros(size(lags))
        ax = OffsetArrays.no_offset_view(axes(x,1))
        lx = length(x)
        for i in 1:lx
            for (ik,k) in enumerate(lags)
                i > lx-k && continue
                if (pred(x[ax[i]])::Bool || pred(x[ax[i+k]])::Bool); cnt[ik] += 1; end
            end
        end
        cnt
    end
    x = repeat([1,2,missing,missing,5], 500)
    lags = 0:4
    @btime count_forlags($ismissing, $x, $lags)
    @btime count_forlags_slower($ismissing, $x, $lags)
    # the version with calling subfunction is faster
    @code_warntype count_forlags(ismissing, x, lags)
    @code_warntype LogNormals.count_forlags2(ismissing, x, lags)
    @code_warntype LogNormals.count_forlag(ismissing, x, 2)
end
