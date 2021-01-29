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

using Lognormals
using StaticArrays, Distributions


D = fit(LogNormal, @qs_muu(3,9))

