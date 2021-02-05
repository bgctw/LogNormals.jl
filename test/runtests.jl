using Test, Documenter, Distributions, LogNormals 

# two moments
m = Moments(1,0.5)
@test n_moments(typeof(m)) == 2
@test n_moments(m) == 2
#@test n_moments(typeof(M)) == 2
@test mean(m) == m[1] == 1
@test var(m) == m[2] == 0.5
@test std(m) == sqrt(0.5)
@test_throws Exception skewness(m)
@test_throws Exception kurtosis(m)

# no moments
m = Moments()
@test n_moments(m) == 0
@test_throws Exception mean(m)
@test_throws Exception m[1]

# moments function of Distribution
d = Normal(2,5)
function testMoments(d)
    m = moments(d, Val(4))
    @test [m[i] for i in 1:n_moments(m)] == [mean(d), var(d), skewness(d), kurtosis(d)]
    m = moments(d, Val(3))
    @test [m[i] for i in 1:n_moments(m)] == [mean(d), var(d), skewness(d)]
    m = moments(d, Val(2))
    @test [m[i] for i in 1:n_moments(m)] == [mean(d), var(d)]
    m = moments(d, Val(1))
    @test [m[i] for i in 1:n_moments(m)] == [mean(d)]
    m = moments(d, Val(0))
    @test [m[i] for i in 1:n_moments(m)] == []
end
testMoments(d)
testMoments(LogNormal(2,5))
# mean not defined fro LogitNormal
@test_throws MethodError testMoments(LogitNormal(2,5))


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
dn = fit(Normal, qpl, qpu)
@test quantile.(dn, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]
dn = fit(Normal, qpu, qpl) # sort
@test quantile.(dn, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]

#println("fitting lognormal")
d = fit(LogNormal, qpl, qpu);
@test quantile.(d, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]
d = fit(LogNormal, qpu, qpl) # sort
@test quantile.(d, [qpl.p, qpu.p]) ≈ [qpl.q, qpu.q]

#println("Approximate Normal by lognormal")
d = fit(LogNormal, moments(dn));
@test mean(d) == mean(dn) && var(d) == var(dn)

if (false) # only interactively
    using StatsPlots
    plot(d); 
    plot!(dn, linetype = :line)
    vline!([mean(d)])
end


#println("fit to quantilePoint and mean")
d = LogNormal(1,1)
m = log(mean(d))
#@macroexpand @qp(0.95,quantile(D,0.95))
qp = @qp(0.95,quantile(d,0.95))
dfit = LogNormals.fit_mean_quantile(LogNormal, mean(d), qp)
@test dfit ≈ d
dfit = fit(LogNormal, mean(d), qp, Val(:mean))
@test dfit ≈ d
# with lower quantile
qp = @qp(0.05,quantile(d,0.05))
dfit = LogNormals.fit_mean_quantile(LogNormal, mean(d), qp)
@test dfit ≈ d
# very close to mean can give very different results:
qp = @qp(0.95,mean(d)-1e-4)
dfit = LogNormals.fit_mean_quantile(LogNormal, mean(d), qp)
@test mean(dfit) ≈ mean(d) && quantile(dfit, qp.p) ≈ qp.q
if (false) # only interactively
    using StatsPlots
    plot(d); plot!(dfit, linetype = :line)
    vline!([mean(d)])
end

#println("fit to quantilePoint and mode")
d = LogNormal(1,1)
m = log(mode(d))
#@macroexpand @qp(0.95,quantile(D,0.95))
qp = @qp(0.95,quantile(d,0.95))
dfit = LogNormals.fit_mode_quantile(LogNormal, mode(d), qp)
@test dfit ≈ d
dfit = fit(LogNormal, mode(d), qp, Val(:mode))
@test dfit ≈ d
# with lower quantile
qp = @qp(0.025,quantile(d,0.025))
dfit = LogNormals.fit_mode_quantile(LogNormal, mode(d), qp)
@test mode(dfit) ≈ mode(d) && quantile(dfit, qp.p) ≈ qp.q
if (false) # only interactively
    using StatsPlots
    plot(d); plot!(dfit, linetype = :line, xlim=(0,quantile(d, 0.975)))
    vline!([mode(d), mode(dfit)])
end

#println("fit to quantilePoint and median")
d = LogNormal(1,1)
qp = @qp(0.95,quantile(d,0.95))
dfit = fit(LogNormal, median(d), qp, Val(:median))
@test dfit ≈ d

# Normal
d = Normal(3,2)
dfit = fit(Normal, moments(d))
@test dfit ≈ d
qp = @qp(0.95,quantile(d,0.95))
dfit = fit(Normal, median(d), qp, Val(:median))
@test dfit ≈ d
dfit = fit(Normal, mode(d), qp, Val(:mode))
@test dfit ≈ d

# testing examples
# make sure to not test for error. This does not work in test, because error compromises former output
# better keep docu and tests separated, but test execution of examples too.
#DocMeta.setdocmeta!(LogNormals, :DocTestSetup, :(using Distributions,LogNormals); recursive=true)
#doctest(LogNormals, manual = false)


