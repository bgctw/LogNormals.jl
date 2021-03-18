using LogNormals
using Test, Distributions, StatsBase, Missings, MissingStrategies, Random
using OffsetArrays, RecursiveArrayTools, LinearAlgebra
using DistributionVectors # cormatri_for_acf
using Unitful

@testset "sem_cor" begin
    acf0 = [1,0.4,0.1]
    Sigma = cormatrix_for_acf(100, acf0);
    dmn = MvNormal(ones(100), Symmetric(Sigma));
    Random.seed!(1234)
    a = rand(dmn);    
    am = allowmissing(a); am[3:4] .= missing;
    az = copy(a); az[3:4] .= 0.0;
    amu = am .* u"m"; # unitful missing vector of meters
    @testset "count_forlags" begin
        pred(x_i, x_iplusk) = ismissing(x_i) || ismissing(x_iplusk)
        x = [1,2,missing,missing,5,6]
        lags = 0:5
        nmiss = @inferred count_forlags(pred, x, lags)
        @test nmiss == [2,3,4,2,0,0]
        x = [1,2,missing,missing,5]
        lags = 0:4
        nmiss = @inferred count_forlags(pred, x, lags)
        @test nmiss == [2,3,3,1,0]
        nmissu = @inferred count_forlags(pred, x.*u"m", lags)
        @test nmissu == nmiss
        # test for custom indices
        nmiss2 = count_forlags(pred, OffsetArray(x,-1), lags)
        @test nmiss2 == nmiss
        # views still work and larger lags return 0
        nmiss3 = @inferred count_forlags(pred, view(x,1:2:5), lags)
        @test nmiss3 == [1,2,0,0,0]
    end;
    @testset "autocorr with missings" begin
        # for comparability center all around nomissing mean
        meana = mean(a)
        lags = 0:8
        acfe1 = @inferred autocor(a, lags) # StatsBase still works
        acfe2 = @inferred autocor(a, lags, PassMissing()) # can supply Missing clause
        @test acfe1 == acfe2
        # get only missing if not explicitly requesting SkipMissing or ExactMissing
        ismissing(@inferred(Vector{Float64},autocor(am, lags, PassMissing(); demean=false)))
        ismissing(@inferred(Vector{Float64},autocor(am, lags)))
        ismissing(@inferred(Vector{Float64},autocor(am, lags; demean=false)))
        acfe = @inferred autocor(a.-meana, lags; demean=false);
        acfem = @inferred(autocor(am.-meana, lags, ExactMissing(); demean=false))
        acfes = @inferred autocor((a.-meana)[4:end], lags; demean=false) # different nobs
        azdemean = copy(a.-meana); azdemean[ismissing.(am)] .= 0.0
        acfez = @inferred autocor(azdemean, lags; demean=false);
        acfemf = @inferred autocor(am.-meana, lags, SkipMissing(); demean=false);
        hcat(acfe, acfem, acfez, acfemf)[1:4,:]
        # with SkipMissing() rather than ExactMissing(): missing same as vector with zeros
        @test acfemf == acfez
        @test acfem[1] == acfez[1] == 1
        # with excactmissing: correct var(x) by two missings (100-2) 
        # lag1: three missings
        @test acfem[2] ≈ acfez[2]*98/97
        # lag2: four missings
        @test acfem[3] ≈ acfez[3]*98/96
        # lag3: three missings
        @test acfem[4] ≈ acfez[4]*98/97
        # lag4: two missings as in z
        @test acfem[5:end] ≈ acfez[5:end]
        # TODO after StatsBase.autocor can deal with unitful
        #acfemu = @inferred autocor((am.-meana).*u"m", lags; demean=false);
        #
        # variant without lags
        ismissing(@inferred autocor(am))
        ismissing(@inferred autocor(am, PassMissing()))
        @test acfemf ≈ (@inferred autocor(am.-meana, SkipMissing(); demean=false))[1:9]
        @test acfem ≈ (@inferred(autocor(am.-meana, ExactMissing(); demean=false)))[1:9]
        #SimpleTraits.istrait(IsEltypeSuperOfMissing{typeof(am.-meana)})
    end;
    @testset "autocor_effective" begin
        effa = @inferred autocor_effective(a, acf0)
        @test effa == acf0
        effa = @inferred autocor_effective(a, [acf0..., -0.001])
        @test effa == acf0
    end;
    @testset "effective_n_cor" begin
        neff = @inferred effective_n_cor(a, acf0)
        @test neff < length(a)
        @test neff ≈ 50.30 atol=0.01 # regression test
        # call with MissingStrategy but nonmissing type
        neffs = @inferred effective_n_cor(a, acf0, ExactMissing())
        @test neffs == neff
        #
        @test ismissing(@inferred(Missing,effective_n_cor(am, acf0, PassMissing())))
        @test ismissing(@inferred(Missing,effective_n_cor(am, acf0)))
        neffm = @inferred effective_n_cor(am, acf0, ExactMissing())
        @test neffm < neff
        @test neffm ≈ 49.61 atol=0.01 # regression test
        # n smaller than length(acf):
        neff2 = @inferred effective_n_cor(a[1:2], [acf0..., 0.05, 0.05])
        @test neff2 ≈ 1.43 atol=0.01 # regression test
        # no correlation estimate for lag 3
        neff3 = @inferred effective_n_cor(am[1:4], [acf0..., 0.05, 0.05], ExactMissing())
        @test neff3 ≈ 1.43 atol=0.01 # regression test
        #
        # without specifying acf
        ismissing(@inferred(Float64,effective_n_cor(am)))
        neff = @inferred effective_n_cor(am, ExactMissing())
        @test neff == effective_n_cor(am, autocor(am, ExactMissing()), ExactMissing())
    end;
    @testset "var_cor" begin
        va1 = @inferred var_cor(a, acf0, ExactMissing())
        va2 = @inferred var_cor(a, acf0, PassMissing())
        va = @inferred var_cor(a, acf0)
        @test va == va1 == va2
        # need to specify handling of missing otherwise get missing back
        ismissing(@inferred(Float64,var_cor(am, acf0)))
        vam = @inferred var_cor(am, acf0, ExactMissing())
        #@code_warntype var_cor(am, acf0, ms=ExactMissing(); neff=nothing)
        @test isnan(var_cor(am[[3]],acf0, ExactMissing()))
        @test isnan(var_cor(am[2:3],acf0, ExactMissing())) # only one non-missing
        vamf = @inferred var_cor(am, acf0, SkipMissing())
        # larger uncertainty with missings
        @test vam > va
        # not correcting neff for missings overestimates neff
        # and hence underestimates uncertainty
        @test vamf < vam
    end;
    @testset "semcor with same effective acf" begin
        se_a = @inferred sem_cor(a, acf0)
        ismissing(@inferred(Float64, sem_cor(am, acf0)))
        @code_warntype sem_cor(am, acf0, PassMissing())
        @edit sem_cor(am, acf0, PassMissing())
        se_am = @inferred sem_cor(am, acf0, ExactMissing())
        @test isnan(sem_cor(am[[3]],acf0, ExactMissing()))
        se_az = @inferred sem_cor(az, acf0)
        ar = filter(x -> !ismissing(x), am)
        se_ar = @inferred Missing sem_cor(ar, acf0)
        [se_a, se_am, se_az, se_ar]
        # the few missing do not lead to large deviations
        @test se_am ≈ se_a atol=0.02
        # skipping missing values has larger stderror, because shorter series
        @test se_ar > se_a
        # caring for missing pairs in correlation also has larger stderror
        # because some of larger relative error in sum_normals
        @test se_am > se_a 
        # the SkipMissing() case may lead to lower estimates
        # of unceertainty than without missings?
    end;
    @testset "semcor calling without acf" begin
        @test @inferred sem_cor(a, ExactMissing()) ≈ sem_cor(a, autocor_effective(a))
        @test @inferred sem_cor(a, SkipMissing()) ≈ sem_cor(a, autocor_effective(a))
        @test @inferred sem_cor(a, PassMissing()) == sem_cor(a, autocor_effective(a))
        @test @inferred sem_cor(am, ExactMissing()) == sem_cor(am, 
            autocor_effective(am, ExactMissing()), ExactMissing())
        @test @inferred sem_cor(am, SkipMissing()) == sem_cor(am, 
            autocor_effective(am, SkipMissing()), SkipMissing())
    end;
    @testset "semcor with many random missings" begin
        b = allowmissing(a);
        Random.seed!(1234)
        b[sample(axes(a,1),60, replace=false)] .= missing;
        se_a = @inferred sem_cor(a, acf0)
        se_b = @inferred Missing sem_cor(b, acf0, ExactMissing())
        se_be = @inferred Missing sem_cor(b, ExactMissing())
        # uncertainty increases due to missings (fewer obs)
        @test se_b > se_a
        #@test se_be > se_a
        # regression tests with random seed 1234
        @test se_b ≈ 0.18 atol = 0.01
    end;
   @testset "semcor short series" begin
        @test isnan(@inferred(var_cor([1.0])))
        @test isnan(@inferred(sem_cor([1.0])))
        @test @inferred sem_cor(1.0:2) == 0.5
        @test @inferred sem_cor(1.0:3) == sqrt(var(1:3)/3)
        @test @inferred var_cor(fill(1.0,4)) == 0.0 
        @test @inferred sem_cor(fill(1.0,4)) == 0.0 # NA correlation but 0 variance
        @test ismissing(@inferred(Float64, sem_cor([1.0,2,missing])))
        @test ismissing(@inferred(Float64, var_cor([1.0,2,missing])))
        @test @inferred var_cor([1.0,2,missing], ExactMissing()) == 0.5
        @test @inferred sem_cor([1.0,2,missing], ExactMissing()) == sqrt(0.5/2)
        @test isnan(@inferred sem_cor([1,missing,missing], ExactMissing()))
   end;
end;

@testset "autocorr with many missings" begin
    # ovrriding a from sd testset
    nsum = 1_000
    acf0 = [1,0.4,0.1]
    Sigma = cormatrix_for_acf(nsum, acf0);
    dmn = MvNormal(ones(nsum), Symmetric(Sigma));
    a = rand(dmn);
    probgap = 0.4
    clustersize = 4
    # need to average over some clusters to get reliable estimate
    nboot = 30
    resboot = map(1:nboot) do iboot
        igap0 = sample(0:(nsum-clustersize),trunc(Int,nsum*probgap/clustersize), replace=true);
        b = allowmissing(a);
        # outer product see https://stackoverflow.com/a/44592419
        #igap0.+Base.OneTo(clustersize)' 
        b[igap0.+Base.OneTo(clustersize)'] .= missing;
        ae = autocor(b, 1:2, ExactMissing())
        af = autocor(b, 1:2, SkipMissing())
        (ae, af)
    end;
    aeb, afb = VectorOfArray.(vectuptotupvec(resboot));
    ae = mean(aeb; dims=2)
    af = mean(afb; dims=2)
    aorig = autocor(a, 1:2)
    @test ae ≈ aorig atol=0.02
    @test all((aorig .- af)./aorig .> 0.05) # >5% low bias
end;

