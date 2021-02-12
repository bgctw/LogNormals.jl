using LogNormals
using Missings, Distributions, Test

@testset "DistributionVector" begin

#a = similar([(0.0, 0.0)], 3)
mu = [1.1,1.2,1.3]
sigma = [1.01, 1.02, 1.03]
a = collect(zip(mu, sigma))

dv = dv0 = DistributionVector(LogNormal, allowmissing(a))

am = allowmissing(a); am[1] = missing
dvm = DistributionVector(LogNormal, am)

@testset "checking constructors" begin
    # not using allowmissing
    @test_throws ErrorException DistributionVector(LogNormal, a) 
    # not allowing others than tuples
    @test_throws ErrorException DistributionVector(LogNormal, allowmissing([1,2])) 
    # not allowing for Tuples of different length
    # already coverey by Vector{T}
end;    

@testset "iterator nonmissing" begin
    @test @inferred length(dv) == 3
    d = @inferred Missing dv[1]
    @test isa(d, LogNormal{Float64})
    @test params(d) == a[1]
    darr = [d for d in dv]
end;    

@testset "iterator missing" begin
    @test @inferred length(dvm) == 3
    ismissing(@inferred LogNormal dvm[1])
    d = @inferred Missing dvm[2]
    @test isa(d, LogNormal{Float64})
    @test params(d) == am[2]
    darr = [d for d in dvm]
    ismissing(darr[1])
end; 

@testset "constructor with several Distributions" begin
    d1 = LogNormal(log(110), 0.25)
    d2 = LogNormal(log(100), 0.15)
    dv = DistributionVector(d1, d2);
    @test @inferred Missing dv[1] == d1
    # empty not allowed
    @test_throws Exception DistributionVector()
    # different types not allowed
    @test_throws MethodError DistributionVector(d1, d2, Normal());
    # with missing 
    dv = DistributionVector(d1, d2, missing);
    @test ismissing(dv[3])
end;

@testset "constructor with parameter vectors" begin
    dv = DistributionVector(LogNormal, mu, sigma)
    d1 = dv[1]
    @test params(d1) == (mu[1], sigma[1])
    # when one parameter has missing, the entire tuple must be set to missing
    mum = allowmissing(mu); mum[1] = missing
    dv = DistributionVector(LogNormal, mum, sigma)
    @test ismissing(dv[1])
end;

@testset "accessing parameters as array" begin
    @test params(dv0,1) == mu
    @test params(dv0,2) == sigma
    # missing
    @test ismissing(params(dvm, 1)[1]) 
    @test params(dvm, 1)[2:3] == mu[2:3]
    @test_throws BoundsError params(dv0,3)
end;

@testset "Biomial distribution with different parameter types" begin
    db1 = Binomial(1, 0.5)
    db2 = Binomial(2, 0.4)
    dv = DistributionVector(db1, db2, missing)
    params(dv, 1)
end;

@testset "Multivaraite distribution with complex parameter types" begin
    #dmn1 = MvNormal(3, 1) # does not work because its of different specific type
    dmn1 = MvNormal([0,0,0], 1)
    dmn2 = MvNormal([1,1,1], 2)
    #params(dmn1), params(dmn2)
    dv = DistributionVector(dmn1, dmn2, missing)
    @test nonmissingtype(eltype(params(dv, 1))) <: AbstractVector
    @test nonmissingtype(eltype(params(dv, 2))) <: AbstractMatrix
end;

end; # testset "DistributionVector"