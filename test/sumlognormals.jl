using LogNormals
using Test, Distributions, LinearAlgebra, Missings

# @testset "two vars uncorrelated" begin
#     # generate nSample values of two lognormal random variables
#     d1 = LogNormal(log(110), 0.25)
#     d2 = LogNormal(log(100), 0.15)
#     dv = DistributionVector(d1, d2);
#     dsum = sum(dv)
#     @testset "no missings" begin
#         @test params(dsum)[1] ≈ log(210) rtol = 0.02 
#         # regression to previous result checked with random numbers below
#         @test exp(params(dsum)[2]) ≈ 1.16087 rtol = 0.02
#     end;
#     @testset "with missings" begin
#         dv = DistributionVector(d1, missing);
#         @test coalesce.(collect(dv), LogNormal()) == [d1, LogNormal()]
#         @test_throws Exception dsum2 = sum(dv) # without skipmissings
#         dsum2 = sum(skipmissing(dv))
#         @test dsum2 == d1
#         dsum3 = sum(dv; skipmissings = Val(true))
#         @test dsum3 == d1
#     end
# end;

# @testset "few correlated vars" begin
#   mu = log.([110,100,80,120,160.0])
#   sigma = log.([1.2,1.5,1.1,1.3,1.1])
#   acf1 = [0.4,0.1]
#   n = length(mu)
#   corrM = cormatrix_for_acf(n, acf1)
#   #
#   @testset "matrix without missing" begin
#     dv = DistributionVector(LogNormal, mu, sigma)
#     dsum = sum(dv, corrM )
#     # regression test TODO check with random numbers
#     @test σstar(dsum) ≈ 1.133632 rtol = 0.02
#     # acf variant
#     dsum3 = sum(dv, acf1)
#     @test all(params(dsum3) .≈ params(dsum))
#   end
#   @testset "matrix with missing" begin
#     mum = allowmissing(mu); mum[1] = missing
#     dsm = DistributionVector(LogNormal, mum, sigma)
#     @test_throws ErrorException dsumm = sum(dsm, corrM )
#     dsumm = sum(dsm, corrM; skipmissings = Val(true) )
#     # regression test TODO
#     @test σstar(dsumm) ≈ 1.15 rtol = 0.02
#     @test_throws ErrorException dsumm2 = sum(dsm, acf1)
#     dsumm2 = sum(dsm, acf1; skipmissings = Val(true))
#     @test all(params(dsumm2) .≈ params(dsumm))
#   #
#   # TODO check by random numbers
#   #nSample = 500
#   #Sigma = Diagonal(sigma) * corrM * Diagonal(sigma)
#   #rM = rand(MvNormal(mu, Symmetric(Sigma)), nSample)
# end;  

# function benchmarkSums()
#     using BenchmarkTools
#     nrep = 30
#     mu = log.(rand(Normal(120, 10), nrep));
#     sigma = log.(rand(Normal(1.2, 0.1), nrep));
#     acf1 = vcat([0.8,0.3,0.1], repeat([0.05], 20));
#     dv = DistributionVector(LogNormal, mu, sigma)
#     @btime dsum_v = sum($dv, $acf1, skipmissings = Val(true), method = Val(:vector) )
#     @btime dsum_m = sum($dv, $acf1, skipmissings = Val(true), method = Val(:bandedmatrix) )
#     # the Bandematrix based is faster
#     #
#     # check type system
#     storage = allowmissing(similar(mu));
#     corMa = cormatrix_for_acf(length(dv),acf1);
#     @code_warntype sum_lognormals!(storage, dv, acf1, skipmissings = Val(true) )
#     @code_warntype sum(dv, acf1, skipmissings = Val(true), method = Val(:vector) )
#     @code_warntype sum_lognormals!(storage, dv, corMa, skipmissings = Val(true) )
#     @code_warntype sum(dv, acf1, skipmissings = Val(true), method = Val(:bandedmatrix) )
#     @btime sum_lognormals!($storage, $dv, $corMa, skipmissings = Val(true) )
#     #
#     #
#     #
#     # repeat with missings
#     mum = allowmissing(copy(mu)); mum[1] = missing
#     dv = DistributionVector(LogNormal, mum, sigma)
#     @btime dsum_v = sum($dv, $acf1, skipmissings = Val(true), method = Val(:vector) )
#     @btime dsum_m = sum($dv, $acf1, skipmissings = Val(true), method = Val(:bandedmatrix) )

#     # try allocating instead of view (replace line of view_nonmissing)
#     # allocating is faster
#     function sum_lognormals2!(S, dv, corr::AbstractMatrix; 
#         skipmissings::Val{l} = Val(false)) where l
#         parms = params(dv)
#         μ = @view parms[1,:]
#         σ = @view parms[2,:]
#         # S = allowmissing(similar(μ))
#         @. S = exp(μ + abs2(σ)/2)
#         nmissing = count(ismissing, S)
#         anymissing = nmissing != 0
#         skipmissings != Val(true) && anymissing && error(
#              "Found missing values. Use argument 'skipmissings = Val(true)' to sum over nonmissing.")
#         Ssum::nonmissingtype(eltype(S)) = sum(skipmissing(S))
#         @. S = σ * S  # do only after Ssum
#         # setting S to zero results in summing zero for missing records
#         # which is the same as filtering both S and corr
#         anymissing && replace!(S, missing => 0.0)
#         #s = transpose(disallowmissing(S)) * corr * disallowmissing(S)
#         #Spure = view_nonmissing(S) # non-allocating
#         Spure = disallowmissing(S) # allocating
#         s = transpose(Spure) * corr * Spure
#         σ2eff = s/abs2(Ssum)
#         μ_sum = log(Ssum) - σ2eff/2
#         #@show Ssum, s, length(S) - nmissing
#         LogNormal(μ_sum, √σ2eff)  
#     end
#     @btime sum_lognormals!($storage, $dv, $corMa, skipmissings = Val(true) )
#     @btime sum_lognormals2!($storage, $dv, $corMa, skipmissings = Val(true) )

    
    
# end


# # function bootstrap_sums_lognormal()
# #     nObs = 5
# #     xTrue = fill(10.0, nObs)
# #     sigmaStar = fill(1.5, nObs) # multiplicative stddev of 1.5
# #     dv = fit.(LogNormal, xTrue, Σstar.(sigmaStar))
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
# #       xObsN[:,i] = rand.(dv,)
# #     end
# #     sums = sum(xObsN, dims = 1)
# #     mean(sums), std(sums)
# #     coefSum = sum_lognormals.( (params.(dv))... )
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
  

