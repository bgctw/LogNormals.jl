plusTwo(3)

pwdDirs = split(pwd(),"/")
dup2 = pwdDirs[end-2:end-1]
if ((dup2 == ["julia","dev"]) | (dup2 == ["storage","julia"]) | (dup2 == ["twutz","julia"]))
  println("Activating current directory. ")
  using Pkg;  Pkg.activate(".")
  println("using Revise")
  using Revise
end

M = Moments(3)

a = [1,2]

using StaticArrays, Distributions
test(a::T...) where T <: Number= SVector{length(a)}(a...)
sa = test(1,3)

abstract type AbstractMoments{N} end
Base.length(::Type{AbstractMoments{N}}) where N = N
Base.length(M::AbstractMoments{N}) where N = N
Distributions.mean(M::AbstractMoments) = length(M) >= 1 ? M[1] : missing
Distributions.var(M::AbstractMoments) = length(M) >= 2 ? M[2] : missing
Distributions.skewness(M::AbstractMoments) = length(M) >= 3 ? M[3] : missing
Distributions.kurtosis(M::AbstractMoments) = length(M) >= 4 ? M[4] : missing

struct Moments{N,T} <: AbstractMoments{N}
    all::SVector{N,T}
end
Moments(x::T...) where T<:Number = Moments(SVector{length(x)}(x...))
Base.getindex(M::Moments, i) = length(M) >= i ? M.all[i] : missing
Base.convert(::Type{AbstractArray}, M::Moments) = M.all

Moments() = Moments(SA[])

Moments(x::T...) where T<:Number = Moments(SVector{length(x)}(x...))
Moments(sa)
Moments(3,1)
Moments(3)