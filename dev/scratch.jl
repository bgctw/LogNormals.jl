"""
    plusTwo(x)

Sum the numeric "2" to whatever it receives as input

A more detailed explanation can go here, although I guess it is not needed in this case

# Arguments
* `x`: The amount to which we want to add 2

# Notes
* Notes can go here

# Examples
```julia
julia> five = plusTwo(3)
5
```
"""
plusTwo(x) = return x+2



M = Moments(3)

a = [1,2]

using LogNormals
using StaticArrays, Distributions


D = fit(LogNormal, @qp_m(3), @qp_uu(9))

function StatsBase.fit(::Type{LogNormal}; mean::Real, mult_std::Real)
    sigma = log(mult_std)
    mu = log(mean) - sigma*sigma/2
    LogNormal(mu, sigma)
end

mult_std(D::LogNormal) = exp(params(D)[2])

@test mean(D) ≈ 4 && mult_std(D) ≈ exp(1)


function fitstats_docu()
    # to specified moments
    d = fit(LogNormal, Moments(3.0,4.0))
    (mean(d), var(d)) .≈ (3.0, 4.0)

    # to mean and upper quantile point
    d = fit(LogNormal, 3, @qp_uu(8))
    (mean(d), quantile(d, 0.975)) .≈ (3.0, 8.0)
    
    # to mode and upper quantile point
    d = fit(LogNormal, 3, @qp_uu(8), Val(:mode))
    (mode(d), quantile(d, 0.975)) .≈ (3.0, 8.0)

    # to two quantiles, i.e confidence range
    d = fit(LogNormal, @qp_ll(1.0), @qp_uu(8))
    (quantile(d, 0.025), quantile(d, 0.975)) .≈ (1.0, 8.0)

    # approximate a different distribution by matching moments
    dn = Normal(3,2)
    d = fit(LogNormal, moments(dn))
    (mean(d), var(d)) .≈ (3.0, 4.0)
end

# some occurence of true
# https://stackoverflow.com/questions/406230/regular-expression-to-match-a-line-that-doesnt-contain-a-word
occursin.(r"^((?!true|\(true).)*$",["true", "(true, true)", "not true", "other", ""])
# true or (true at the beginning
occursin.(r"^(?!true|\(true).",["true", "(true, true)", "not true", "other", ""])


function experimentAR1()
    using StatsBase
    # generate autocorrelated time series
    using ARFIMA
    N, σ = 1_000, 0.5
    x = arfima(N, σ, nothing, SVector(0.6))
    # re-estimate autocorrelation
    autocor(x, 1:10)
end

tvec = allowmissing([(rand(),rand()) for i=1:6]);
tvec[1] = missing
c = mappedarray(x-> ismissing(x) ? missing : first(x),tvec);
c[1]
c[1:2]
c[[1,2]]

using BenchmarkTools: @btime
using Test


x = [missing, 1.0];
y = [missing, 2.0];
f(x,y) = x[1] * y[1] * x[2] * y[2];
@btime f($x,$y)

function f1(x,y)
    s = zero(eltype(x))
    for i in axes(x,1), j in axes(y,1)
        sij = x[i] * x[j] * y[i] * y[j]
        if !ismissing(sij) 
            s += sij
        end
    end
    s
end
@btime f1($x,$y)

f2(x,y) = sum(skipmissing(x[i] * x[j] + y[i] * y[j] for i in axes(x,1), j in axes(y,1)))
@btime f2($x,$y)


mul4(x1,x2,x3,x4) = x1 * x2 * x3 * x4
function f(x,y)
    s = zero(nonmissingtype(eltype(x)))
    for i in axes(x,1), j in axes(y,1)
        sij = passmissing(+)(x[i] * x[j],  y[i] * y[j])
        if !ismissing(sij) 
            s += sij::nonmissingtype(eltype(x))
        end
    end
    s
end
@code_warntype f(x,y)
@time f(x,y); @time f(x,y)

function f(x,y)
    T = nonmissingtype(eltype(x))
    s = zero(T)
    for i in axes(x,1), j in axes(y,1)
        #if !any(ismissing.((x[i], x[j], y[i], y[j])))
        if !ismissing(x[i]) && !ismissing(x[j]) && !ismissing(y[i]) && !ismissing(y[i]) 
            sij = x[i] * x[j] + y[i] * y[j]
            s += sij::T
        end
    end
    s
end
#@code_warntype f(x,y)
@time f(x,y); @time f(x,y)


f2(x,y) = sum(skipmissing(x[i] * x[j] * y[i] * y[j] for i in axes(x,1), j in axes(y,1)))
@btime f2($x,$y)



f() = *(x[1], y[1], x[2], y[2])
@time f(); @time f()

# fast compiler code works on value-type keyword argument
function f(;sk::Val{B} = Val(true)) where B
    B==true && return(1)
    2
end
f()
@code_llvm f()
@code_llvm f(sk=Val(false))

#modifying several arrays with one filter in the same loop
function f!(a,b,F)
    for i in Iterators.filter(F, 1:length(a))
        a[i] = 0
        b[i] = 1
    end
    a,b
end
a = collect(1:5)
b = collect(1:5)
isgapfilled = falses(5); isgapfilled[3] = true
f!(a,b, (i -> isgapfilled[i]))


function f(dv::AbstractDistributionVector{D}) where {D<:Distribution} 
    x1 = rand(first(skipmissing(dv)))
    fmiss(x) = ismissing(x) ? missing : rand(x)
    allowmissing(fmiss.(dv))::Vector{Union{Missing,typeof(x1)}} 
end
f(dv)
@code_warntype f(dv)


x = [missing, 1.0];
y = [missing, 2.0]
g(x) = 2 * x
function f(x)
    #collect(eltype(x), g.(x))::Vector{eltype(x)}
    typeof(x)(g.(x))::typeof(x)
end
f(x)
@inferred f(x)

using StaticArrays
xs = SVector(missing , 2.0)
@inferred f(xs)
