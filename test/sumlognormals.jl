using LogNormals
using Test, Distributions, LinearAlgebra, BandedMatrices, Missings

if false
end

@testset "DistributionSequence" begin
    "constructor with type and matrix" 
    mu = [1.1,1.2,1.3]
    sigma = [1.01, 1.02, 1.03]
    # makes a copy to support missing
    ds = ds1 = DistributionSequence(LogNormal, transpose(hcat(mu, sigma)))
    # must be of type Distribution
    @test_throws Exception DistributionSequence(String, hcat(mu, sigma))
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
        bn = [1,1]
        bp = [0.5, 0.7]
        dsb = DistributionSequence(LogNormal, bn, bp)
    end
    @testset "missings" begin
        mu = [1.1,1.2,missing]
        sigma = [1.01, 1.02, 1.03]
        dsm = DistributionSequence(LogNormal, transpose(hcat(mu, sigma)))
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
        @test_throws MethodError ds2 = DistributionSequence(db1, db2) #nparams not defined
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
  corrM = BandedMatrix{Float64}(undef, (n,n), (2,2))
  corrM[band(0)] .= 1
  for i in 1:length(acf1)
    corrM[band(i)] .= corrM[band(-i)] .= acf1[i]
  end
  corrM = Symmetric(corrM)
  det(Array(corrM))
  nSample = 500
  #Sigma = Diagonal(abs2.(sigma)) * corrM
  Sigma = Diagonal(sigma) * corrM * Diagonal(sigma)
  #
  ds = DistributionSequence(LogNormal, mu, sigma)
  dsum = sum(ds, corrM )
  # regression test TODO check with random numbers
  @test σstar(dsum) ≈ 1.133632 rtol = 0.02
  #
  # repeat with explicitely constraining correlation length
  dsum2 = sum(ds, corrM, corrlength = length(acf1) )
  @test dsum2 == dsum
  #
  # missings
  mum = allowmissing(mu); mum[1] = missing
  dsm = DistributionSequence(LogNormal, mum, sigma)
  @test_throws ErrorException dsumm = sum(dsm, corrM )
  dsumm = sum(dsm, corrM; skipmissings = Val(true) )
  # regression test
  @test σstar(dsumm) ≈ 1.05 rtol = 0.02
  #
  # TODO check by random numbers
  #rM = rand(MvNormal(mu, Symmetric(Sigma)), nSample)
end;  


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
  

