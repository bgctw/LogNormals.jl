using Test, Distributions, Lognormals 

out = plusTwo(3)
@test out == 5

# two moments
M = Moments(1,0.5)
@test length(typeof(M)) == 2
@test length(M) == 2
#@test length(typeof(M)) == 2
@test mean(M) == M[1] == 1
@test var(M) == M[2] == 0.5
@test std(M) == sqrt(0.5)
@test_throws Exception skewness(M)
@test_throws Exception kurtosis(M)
typeof(convert(AbstractArray,M)) <: AbstractArray

# no moments
M = Moments()
@test length(M) == 0
@test_throws Exception mean(M)
@test_throws Exception M[1]
typeof(convert(AbstractArray,M)) <: AbstractArray

# distribution parameters from moments
D = LogNormal(1,0.6)
M = Moments(mean(D), var(D))
Dfit = fit(LogNormal, M)
@test D ≈ Dfit

# handle not giving variance
@test_throws Exception D = fit(LogNormal, Moments(3.2))