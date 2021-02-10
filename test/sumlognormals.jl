using Test, Distributions, LogNormals, LinearAlgebra, BandedMatrices, Missings

if false
end

@testset "two vars coef" begin
  # generate nSample values of two lognormal random variables
  mu1 = log(110)
  mu2 = log(100)
  sigma1 = 0.25
  sigma2 = 0.15
  #estimateSumLognormalBenchmark( c(mu1,mu2), c(sigma1,sigma2) )
  dsum = sum_lognormals( [mu1,mu2], [sigma1,sigma2] )
  @test params(dsum)[1] ≈ log(210) rtol = 0.02 
  # regression to previous result checked with random numbers belowq
  @test exp(params(dsum)[2]) ≈ 1.16087 rtol = 0.02
end;

@testset "two vars" begin
    # generate nSample values of two lognormal random variables
    d1 = LogNormal(log(110), 0.25)
    d2 = LogNormal(log(100), 0.15)
    dsum = sum_lognormals( [d1, d2] )
    @testset "no missings" begin
        @test params(dsum)[1] ≈ log(210) rtol = 0.02 
        # regression to previous result checked with random numbers belowq
        @test exp(params(dsum)[2]) ≈ 1.16087 rtol = 0.02
    end;
    @testset "with missings" begin
        ds = [d1, d2, missing]
        length_itr(skipmissing(ds))
        dsum2 = sum_lognormals(skipmissing(ds))
        @test dsum2 == dsum
    end
end;

@testset "few Vars" begin
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
  rM = rand(MvNormal(mu, Symmetric(Sigma)), nSample)
  #
  dsum = sum_lognormals(mu, sigma, corrM )
  @test exp(params(dsum)[2]) ≈ 1.133632 rtol = 0.02
  #
  # repeat with explicitely constraining correlation length
  dsum = sum_lognormals(mu, sigma, corrM; corrlength = length(acf1) )
  @test exp(params(dsum)[2]) ≈ 1.133632 rtol = 0.02
end;  

@testset "handle missings uncorrelated" begin
  mu = log.([110,100,80.0])
  sigma = log.(repeat([1.2], length(mu)))
  # one finite case
  sigma1 = allowmissing(sigma); sigma1[2:end] .= missing;
  @test_throws Exception sum_lognormals(mu, sigma1)
  dsum = sum_lognormals(mu, sigma1; skipmissings = Val(true))
  @test all(params(dsum) .≈ (mu[1], sigma[1]))
  # no finite case
  sigma0 = allowmissing(sigma); sigma0[:] .= missing;
  @test_throws Exception coefSum = sum_lognormals(mu, sigma0)
  @test_throws Exception sum_lognormals(mu, sigma0; skipmissings = Val(true))
end;

@testset "handle missings correlated" begin
    mu = log.([110,100,80.0])
    sigma = log.(repeat([1.2], length(mu)))
    corrM = [1   0.4  0
             0.4 1    0.4
             0   0.4  1]
    # one finite case
    sigma1 = allowmissing(sigma); sigma1[2:end] .= missing;
    @test_throws Exception sum_lognormals(mu, sigma1, corrM)
    dsum = sum_lognormals(mu, sigma1, corrM; skipmissings = Val(true))
    @test all(params(dsum) .≈ (mu[1], sigma[1]))
    # no finite case
    sigma0 = allowmissing(sigma); sigma0[:] .= missing;
    @test_throws Exception coefSum = sum_lognormals(mu, sigma0)
    @test_throws ErrorException sum_lognormals(mu, sigma0; skipmissings = Val(true))
end;

# function bootstrap_sums_lognormal()
#     nObs = 5
#     xTrue = fill(10.0, nObs)
#     sigmaStar = fill(1.5, nObs) # multiplicative stddev of 1.5
#     ds = fit.(LogNormal, xTrue, Σstar.(sigmaStar))
#     # generate observations with correlated errors
#     # # acf1 <- c(0.4,0.1)
#     # # corrM <- setMatrixOffDiagonals(
#     # #   diag(nrow = nObs), value = acf1, isSymmetric = TRUE)
#     # xObsN <- exp(mvtnorm::rmvnorm(
#     #   100, mean = theta[,1]
#     #   , sigma = diag(theta[,2]) %*% corrM %*% diag(theta[,2])))
#     nRep = 30
#     #nRep = 1e6
#     xObsN = zeros(nObs,nRep)
#     for i in 1:nRep
#       xObsN[:,i] = rand.(ds,)
#     end
#     sums = sum(xObsN, dims = 1)
#     mean(sums), std(sums)
#     coefSum = sum_lognormals.( (params.(ds))... )
#     moments = getLognormMoments(coefSum[["mu"]], coefSum[["sigma"]])
#     c(sqrt(as.vector(moments[,"var"])), sd(sums))
#     expect_equal(sqrt(as.vector(moments[,"var"])), sd(sums), tolerance = 0.5)
#     expect_true((moments[,"mean"] - sum(xTrue)) < 4*sqrt(moments[,"var"]))
#     expect_equal(coefSum, coefSample, tolerance = 0.01)
#     .tmp.f <- function(){
#       plot(ecdf(sums))
#       q = seq(qlnorm(0.025, coefSum[1], coefSum[2]), qlnorm(0.975, coefSum[1], coefSum[2]), length.out = 201)
#       cdf = plnorm(q, coefSum[1], coefSum[2])
#       lines(cdf ~ q, col = "blue")
#       #
#       plot(density(sums))
#       lines(dlnorm(q, coefSum[1], coefSum[2]) ~ q, col = "blue")
#     }
  
# end
  