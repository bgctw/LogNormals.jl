using Test, Distributions, Lognormals # This load both the test suite and our MyAwesomePackage

out = plusTwo(3)
@test out == 5               # This is the actual test condition. You can add as many tests as you wish.

# two moments
M = Moments(1,0.5)
@test length(M) == 2
#@test length(typeof(M)) == 2
@test mean(M) == M[1] == 1
@test var(M) == M[2] == 0.5
@test ismissing(skewness(M)) & ismissing(M[3])
typeof(convert(AbstractArray,M)) <: AbstractArray

# no moments
M = Moments()
@test length(M) == 0
@test ismissing(mean(M)) & ismissing(M[1])
typeof(convert(AbstractArray,M)) <: AbstractArray

