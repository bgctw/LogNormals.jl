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

using ResultTypes
function integer_division(x::Int, y::Int)::Result{Int, DivideError}
    if y == 0
        return DivideError()
    else
        return div(x, y)
    end
end

x = integer_division(3,2)
x == 1