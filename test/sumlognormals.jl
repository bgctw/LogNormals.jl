using LogNormals
using Test, Distributions, LinearAlgebra, BandedMatrices, Missings

if false
end

@testset "DistributionSequence" begin
    "constructor with type and matrix" 
    mu = [1.1,1.2,1.3]
    sigma = [1.01, 1.02, 1.03]
    # makes a copy to support missing
    ds = ds1 = DistributionSequence(LogNormal, allowmissing(transpose(hcat(mu, sigma))))
    # must be of type Distribution
    @test_throws Exception DistributionSequence(String, hcat(mu, sigma))
    # must allow missings
    @test_throws Exception DistributionSequence(LogNormal, transpose(hcat(mu, sigma)))
    @test eltype(ds) == Union{Missing,LogNormal}
    #params(ds)
    #params(ds)[:,2]
    @test ds[2] == LogNormal(1.2, 1.02)
    dsa = collect(ds)
    @test @inferred length(dsa) == 3
    @test dsa[2] == LogNormal(1.2, 1.02)
    dsrv = collect(Iterators.reverse(ds))
    @test dsrv[3] == LogNormal(1.1, 1.01)
    #@test @inferred params(ds,2) == [1.2, 1.02]
    @testset "constructor with several Distributions" begin
        d1 = LogNormal(log(110), 0.25)
        d2 = LogNormal(log(100), 0.15)
        ds = DistributionSequence(d1, d2);
        @test @inferred Missing ds[1] == d1
        # empty not allowed
        @test_throws Exception DistributionSequence()
        # different types not allowed
        @test_throws MethodError DistributionSequence(d1, d2, Normal());
        # with missing 
        ds = DistributionSequence(d1, d2, missing);
        @test ismissing(ds[3])
    end;
    @testset "constructor with parameter vectors" begin
        ds2 = DistributionSequence(LogNormal, mu, sigma)
        @test ds2.params == ds1.params
    end
    @testset "missings" begin
        mu = [1.1,1.2,missing]
        sigma = [1.01, 1.02, 1.03]
        dsm = DistributionSequence(LogNormal, allowmissing(transpose(hcat(mu, sigma))))
        @test eltype(dsm) == Union{Missing,LogNormal}
        # if (false)
        #     @code_warntype(length(dsm))
        #     @code_warntype(dsm[2])
        # end
        @test @inferred length(dsm) == 3
        @test @inferred Missing dsm[2] == LogNormal(1.2, 1.02)
        @test ismissing(dsm[3])
        dsm2 = DistributionSequence(LogNormal, mu, sigma)
        @test coalesce.(dsm2.params,0815) == coalesce.(dsm.params, 0815)
    end
    @testset "explicit specification of Npar" begin
        db1 = Binomial()
        db2 = Binomial(1,0.7)
        #ds = (db1, db2)
        @test_throws MethodError ds2 = DistributionSequence(db1, db2) #nparams not defineda = 
        ds2 = @inferred DistributionSequence(db1, db2, npar=Val(2)) 
    end
end;

@testset "two vars uncorrelated" begin
    # generate nSample values of two lognormal random variables
    d1 = LogNormal(log(110), 0.25)
    d2 = LogNormal(log(100), 0.15)
    ds = DistributionSequence(d1, d2);
    dsum = sum(ds)
    @testset "no missings" begin
        @test params(dsum)[1] ≈ log(210) rtol = 0.02 
        # regression to previous result checked with random numbers below
        @test exp(params(dsum)[2]) ≈ 1.16087 rtol = 0.02
    end;
    @testset "with missings" begin
        ds = DistributionSequence(d1, missing);
        @test coalesce.(collect(ds), LogNormal()) == [d1, LogNormal()]
        @test_throws Exception dsum2 = sum(ds) # without skipmissings
        dsum2 = sum(skipmissing(ds))
        @test dsum2 == d1
        dsum3 = sum(ds; skipmissings = Val(true))
        @test dsum3 == d1
    end
end;

@testset "few correlated vars" begin
  mu = log.([110,100,80,120,160.0])
  sigma = log.([1.2,1.5,1.1,1.3,1.1])
  acf1 = [0.4,0.1]
  n = length(mu)
  corrM = cormatrix_for_acf(n, acf1)
  #
  @testset "matrix without missing" begin
    ds = DistributionSequence(LogNormal, mu, sigma)
    dsum = sum(ds, corrM )
    # regression test TODO check with random numbers
    @test σstar(dsum) ≈ 1.133632 rtol = 0.02
    # acf variant
    dsum3 = sum(ds, acf1)
    @test all(params(dsum3) .≈ params(dsum))
  end
  @testset "matrix with missing" begin
    mum = allowmissing(mu); mum[1] = missing
    dsm = DistributionSequence(LogNormal, mum, sigma)
    @test_throws ErrorException dsumm = sum(dsm, corrM )
    dsumm = sum(dsm, corrM; skipmissings = Val(true) )
    # regression test TODO
    @test σstar(dsumm) ≈ 1.15 rtol = 0.02
    @test_throws ErrorException dsumm2 = sum(dsm, acf1)
    dsumm2 = sum(dsm, acf1; skipmissings = Val(true))
    @test all(params(dsumm2) .≈ params(dsumm))
  #
  # TODO check by random numbers
  #nSample = 500
  #Sigma = Diagonal(sigma) * corrM * Diagonal(sigma)
  #rM = rand(MvNormal(mu, Symmetric(Sigma)), nSample)
end;  

function benchmarkSums()
    using BenchmarkTools
    nrep = 30
    mu = log.(rand(Normal(120, 10), nrep));
    sigma = log.(rand(Normal(1.2, 0.1), nrep));
    acf1 = vcat([0.8,0.3,0.1], repeat([0.05], 20));
    ds = DistributionSequence(LogNormal, mu, sigma)
    @btime dsum_v = sum($ds, $acf1, skipmissings = Val(true), method = Val(:vector) )
    @btime dsum_m = sum($ds, $acf1, skipmissings = Val(true), method = Val(:bandedmatrix) )
    # the Bandematrix based is faster
    #
    # check type system
    storage = allowmissing(similar(mu));
    corMa = cormatrix_for_acf(length(ds),acf1);
    @code_warntype sum_lognormals!(storage, ds, acf1, skipmissings = Val(true) )
    @code_warntype sum(ds, acf1, skipmissings = Val(true), method = Val(:vector) )
    @code_warntype sum_lognormals!(storage, ds, corMa, skipmissings = Val(true) )
    @code_warntype sum(ds, acf1, skipmissings = Val(true), method = Val(:bandedmatrix) )
    @btime sum_lognormals!($storage, $ds, $corMa, skipmissings = Val(true) )
    #
    #
    #
    # repeat with missings
    mum = allowmissing(copy(mu)); mum[1] = missing
    ds = DistributionSequence(LogNormal, mum, sigma)
    @btime dsum_v = sum($ds, $acf1, skipmissings = Val(true), method = Val(:vector) )
    @btime dsum_m = sum($ds, $acf1, skipmissings = Val(true), method = Val(:bandedmatrix) )

    # try allocating instead of view (replace line of view_nonmissing)
    # allocating is faster
    function sum_lognormals2!(S, ds, corr::AbstractMatrix; 
        skipmissings::Val{l} = Val(false)) where l
        parms = params(ds)
        μ = @view parms[1,:]
        σ = @view parms[2,:]
        # S = allowmissing(similar(μ))
        @. S = exp(μ + abs2(σ)/2)
        nmissing = count(ismissing, S)
        anymissing = nmissing != 0
        skipmissings != Val(true) && anymissing && error(
             "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
        Ssum::nonmissingtype(eltype(S)) = sum(skipmissing(S))
        @. S = σ * S  # do only after Ssum
        # setting S to zero results in summing zero for missing records
        # which is the same as filtering both S and corr
        anymissing && replace!(S, missing => 0.0)
        #s = transpose(disallowmissing(S)) * corr * disallowmissing(S)
        #Spure = view_nonmissing(S) # non-allocating
        Spure = disallowmissing(S) # allocating
        s = transpose(Spure) * corr * Spure
        σ2eff = s/abs2(Ssum)
        μ_sum = log(Ssum) - σ2eff/2
        #@show Ssum, s, length(S) - nmissing
        LogNormal(μ_sum, √σ2eff)  
    end
    @btime sum_lognormals!($storage, $ds, $corMa, skipmissings = Val(true) )
    @btime sum_lognormals2!($storage, $ds, $corMa, skipmissings = Val(true) )

    
    
end


# # function bootstrap_sums_lognormal()
# #     nObs = 5
# #     xTrue = fill(10.0, nObs)
# #     sigmaStar = fill(1.5, nObs) # multiplicative stddev of 1.5
# #     ds = fit.(LogNormal, xTrue, Σstar.(sigmaStar))
# #     # generate observations with correlated errors
# #     # # acf1 <- c(0.4,0.1)
# #     # # corrM <- setMatrixOffDiagonals(
# #     # #   diag(nrow = nObs), value = acf1, isSymmetric = TRUE)
# #     # xObsN <- exp(mvtnorm::rmvnorm(
# #     #   100, mean = theta[,1]
# #     #   , sigma = diag(theta[,2]) %*% corrM %*% diag(theta[,2])))
# #     nRep = 30
# #     #nRep = 1e6
# #     xObsN = zeros(nObs,nRep)
# #     for i in 1:nRep
# #       xObsN[:,i] = rand.(ds,)
# #     end
# #     sums = sum(xObsN, dims = 1)
# #     mean(sums), std(sums)
# #     coefSum = sum_lognormals.( (params.(ds))... )
# #     moments = getLognormMoments(coefSum[["mu"]], coefSum[["sigma"]])
# #     c(sqrt(as.vector(moments[,"var"])), sd(sums))
# #     expect_equal(sqrt(as.vector(moments[,"var"])), sd(sums), tolerance = 0.5)
# #     expect_true((moments[,"mean"] - sum(xTrue)) < 4*sqrt(moments[,"var"]))
# #     expect_equal(coefSum, coefSample, tolerance = 0.01)
# #     .tmp.f <- function(){
# #       plot(ecdf(sums))
# #       q = seq(qlnorm(0.025, coefSum[1], coefSum[2]), qlnorm(0.975, coefSum[1], coefSum[2]), length.out = 201)
# #       cdf = plnorm(q, coefSum[1], coefSum[2])
# #       lines(cdf ~ q, col = "blue")
# #       #
# #       plot(density(sums))
# #       lines(dlnorm(q, coefSum[1], coefSum[2]) ~ q, col = "blue")
# #     }
  
# # end
  

