using LogNormals
using Test, Distributions, LinearAlgebra, Missings

function boot_dvsum(dv, nboot = 100_000)
    nboot = 100_000
    x = rand(dv, nboot);
    sums = map(sum, x);
    params(fit(LogNormal,sums))
end

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

@testset "sumlognormals" begin

@testset "two vars uncorrelated" begin
    # generate nSample values of two lognormal random variables
    d1 = LogNormal(log(110), 0.25)
    d2 = LogNormal(log(100), 0.15)
    dv = SimpleDistributionVector(d1, d2);
    dsum = @inferred sum(dv)
    @testset "no missings" begin
        @test params(dsum)[1] ≈ log(210) rtol = 0.02 
        # checked with random numbers 
        #boot_dvsum(dv)
        @test params(dsum)[2] ≈ 0.15 rtol = 0.02
    end;
    @testset "with missings" begin
        dv = SimpleDistributionVector(d1, missing);
        @test coalesce.(collect(dv), LogNormal()) == [d1, LogNormal()]
        @test_throws Exception dsum2 = sum(dv) # without skipmissings
        dsum2 = @inferred sum(skipmissing(dv))
        @test dsum2 == d1
        dsum3 = @inferred sum(dv; skipmissings = Val(true))
        @test dsum3 == d1
        #@btime sum(skipmissing($dv)) # does not allocate
    end
    @testset "with gapfilling flag" begin
      dv = SimpleDistributionVector(d1, d2, d1);
      isgapfilled = [true, false, false]
      dsum4 = @inferred sum(dv)
      dsum = @inferred sum(dv, isgapfilled)
      #@code_warntype sum(dv, isgapfilled)
      @test mean(dsum) == mean(dsum4)
      @test std(dsum) > std(dsum4)
    end
    @testset "with missings and gapfilling flag" begin
      dv = SimpleDistributionVector(d1, d2, d1, missing);
      isgapfilled = [true, true, false, false]
      dsum5 = @inferred sum(dv[1:3], isgapfilled[1:3])
      dsum = @inferred sum(dv, isgapfilled; skipmissings = Val(true))
      #@code_warntype sum(dv, isgapfilled; skipmissings = Val(true))
      @test dsum == dsum5
    end
end;

@testset "few correlated vars" begin
  mu = log.([110,100,80,120,160.0])
  sigma = log.([1.2,1.5,1.1,1.3,1.1])
  acf1 = [0.4,0.1]
  n = length(mu)
  corrM = @inferred cormatrix_for_acf(n, acf1)
  dv = SimpleDistributionVector(LogNormal{eltype(mu)}, mu, sigma);
  mum = allowmissing(mu); mum[1] = missing
  dvm = SimpleDistributionVector(LogNormal{eltype(mu)}, mum, sigma)
  #
  @testset "matrix without missing" begin
    dsum = @inferred sum(dv, corrM )
    # checked with random numbers
    #boot_dvsums_acf(dv, acf1)
    @test params(dsum)[2] ≈ 0.128 rtol = 0.02
    @test mean(dsum) ≈ sum(mean.(dv))
    # acf variant
    dsum3 = @inferred sum(dv, acf1)
    @test all(params(dsum3) .≈ params(dsum))
  end
  @testset "matrix with missing" begin
    @test_throws ErrorException dsumm = sum(dvm, corrM )
    #S = similar(mum);
    dsumm = @inferred sum(dvm, corrM; skipmissings = Val(true) )
    # checked from sample
    # boot_dvsums_acf(dvm[2:end], acf1)
    #@test σstar(dsumm) ≈ 1.15 rtol = 0.02
    @test params(dsumm)[2] ≈ 0.141 rtol = 0.02
    @test_throws ErrorException dsumm2 = sum(dvm, acf1)
    dsumm2 = @inferred sum(dvm, acf1; skipmissings = Val(true))
    @test all(params(dsumm2) .≈ params(dsumm))
  end;
  @testset "with gapfilling flag" begin
    isgapfilled = fill(false, length(dv)); isgapfilled[4:end] .= true
    dsum4 = sum(dv, corrM)
    dsum = @inferred sum(dv, corrM, isgapfilled)
    #@code_warntype sum(dv, isgapfilled)
    @test mean(dsum) == mean(dsum4)
    @test std(dsum) > std(dsum4)
  end
  @testset "with missings and gapfilling flag" begin
    isgapfilled = fill(false, length(dvm)); isgapfilled[4:end] .= true
    dsum4 = sum(dvm, corrM; skipmissings = Val(true))
    dsum4b = @inferred sum(dvm, corrM, isgapfilled; skipmissings = Val(true))
    #@code_warntype sum(dvm, isgapfilled)
    @test mean(dsum4b) == mean(dsum4)
    @test std(dsum4b) > std(dsum4)
    # acf variant 
    dsum5b = @inferred sum(dvm, acf1, isgapfilled; skipmissings = Val(true))
    @test dsum5b == dsum4b
  end;
end;  

end; # testset "sumlognormals"
# see also conde in sumlognormals_benchmark.jl
# could not stay here because of btime macro
