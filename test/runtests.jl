using Test, Distributions, LogNormals 

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

# moments function of Distribution
D = Normal(2,5)
function testMoments(D)
    m = moments(D, Val(4))
    @test [m[i] for i in 1:length(m)] == [mean(D), var(D), skewness(D), kurtosis(D)]
    m = moments(D, Val(3))
    @test [m[i] for i in 1:length(m)] == [mean(D), var(D), skewness(D)]
    m = moments(D, Val(2))
    @test [m[i] for i in 1:length(m)] == [mean(D), var(D)]
    m = moments(D, Val(1))
    @test [m[i] for i in 1:length(m)] == [mean(D)]
    m = moments(D, Val(0))
    @test [m[i] for i in 1:length(m)] == []
end
testMoments(D)
testMoments(LogNormal(2,5)


# distribution parameters from moments
D = LogNormal(1,0.6)
M = Moments(mean(D), var(D))
Dfit = fit(LogNormal, M)
@test D ≈ Dfit

# handle not giving variance
@test_throws Exception fit(LogNormal, Moments(3.2))

# Quantile Point and 
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

# macros QuantilePoint
@test @qp(0.4,0.7) == QuantilePoint(0.4,0.7)
@test @qp_ll(0.7) == QuantilePoint(0.025,0.7)
@test @qp_l(0.7) == QuantilePoint(0.05,0.7)
@test @qp_m(0.7) == QuantilePoint(0.5,0.7)
@test @qp_u(0.7) == QuantilePoint(0.95,0.7)
@test @qp_uu(0.7) == QuantilePoint(0.975,0.7)

@test @qs_cf90(0.2,0.7) == Set([@qp_l(0.2),@qp_u(0.7)])
@test @qs_cf95(0.2,0.7) == Set([@qp_ll(0.2),@qp_uu(0.7)])


#println("fitting normal")
qpl = @qp_m(3)
qpu = @qp_u(5)
DN = fit(Normal, qpl, qpu)
@test quantile.(DN, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]
DN = fit(Normal, qpu, qpl) # sort
@test quantile.(DN, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]

#println("fitting lognormal")
D = fit(LogNormal, qpl, qpu);
@test quantile.(D, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]
D = fit(LogNormal, qpu, qpl) # sort
@test quantile.(D, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]

#println("Approximate Normal by lognormal")
D = fit(LogNormal, moments(DN));
@test mean(D) == mean(DN) && var(D) == var(DN)

if (FALSE) # only interactively
    using StatsPlots
    plot(D); 
    plot!(DN, linetype = :line)
    vline!([mean(D)])
end




