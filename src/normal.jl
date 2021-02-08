function Distributions.fit(::Type{Normal}, m::AbstractMoments)
    n_moments(m) >= 2 || error("Need mean and variance to estimate normal")
    return Normal(mean(m),std(m))
end

function Distributions.fit(::Type{Normal}, lower::QuantilePoint, upper::QuantilePoint)
    # https://www.johndcook.com/quantiles_parameters.pdf
    if (upper < lower) 
        lower,upper = (upper,lower)
    end
    qz1 = quantile(Normal(), lower.p)
    qz2 = quantile(Normal(), upper.p)
    dqz = (qz2 - qz1)
    σ = (upper.q - lower.q)/dqz
    μ = (lower.q*qz2 - upper.q*qz1)/dqz
    Normal(μ,σ)
end

fit_mean_quantile(D::Type{Normal}, mean, qp::QuantilePoint) = 
    fit(D, QuantilePoint(0.5, mean), qp)

fit_mode_quantile(D::Type{Normal}, mode, qp::QuantilePoint) = 
    fit(D, QuantilePoint(0.5, mode), qp)


