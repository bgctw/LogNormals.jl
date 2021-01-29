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
@test_throws Exception fit(LogNormal, Moments(3.2))

# Quantile Point and QuantileSet
qp1 = QuantilePoint(0.25,2);
qp2 = QuantilePoint(0.95,7.5);

#println("isless")
# p and q both increasing/decreasing
@test (QuantilePoint(0.4, 2) < QuantilePoint(0.5, 3)) == true
@test (QuantilePoint(0.6, 2) < QuantilePoint(0.5, 1)) == false
# same q may have different p as same samples can be repeated
@test (QuantilePoint(0.4, 2) < QuantilePoint(0.5, 2)) == false
@test (QuantilePoint(0.6, 2) < QuantilePoint(0.5, 2)) == false
# q decreasing but p increasing or vise-versa
@test_throws ErrorException QuantilePoint(0.6, 1) < QuantilePoint(0.5, 2) 
@test_throws ErrorException QuantilePoint(0.4, 2) < QuantilePoint(0.5, 1)
# for same p equal, also q must be equal
@test_throws ErrorException  QuantilePoint(0.5, 1) < QuantilePoint(0.5, 2) 
@test_throws ErrorException  QuantilePoint(0.5, 2) < QuantilePoint(0.5, 1)

#println("modifying QuantilePoint")
qp3 = QuantilePoint(qp1, q = 3);
@test qp3.q == 3
@test qp3.p == qp1.p == 0.25

#println("QuantileSet")
qset = QuantileSet([qp2,qp1]);
@test first(qset)== qp1 # reordered

#println("fitting normal")
DN = Lognormals.normal_from_two_quantiles(qp1,qp2);
@test quantile.(DN, [qp1.p, qp2.p]) ≈ [qp.q for qp in qset]

#println("fitting lognormal")
D = fit(LogNormal, qset);
@test quantile.(D, [qp1.p, qp2.p]) ≈ [qp.q for qp in qset]

#println("fitting end")



