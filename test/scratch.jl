

M = Moments(3)

a = [1,2]

using LogNormals
using StaticArrays, Distributions


D = fit(LogNormal, @qp_m(3), @qp_uu(9))

