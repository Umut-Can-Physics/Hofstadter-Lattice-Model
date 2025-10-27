include("Laughlin Scripts/JacobiThetaFunction.jl")

using Plots

a, b = 0.3, 0.2
τ = 0.5im
N = 50
xs = range(-1, 1, length=200)
ys = range(-1, 1, length=200)

vals = [abs(v(a,b,x + im*y, τ, N)) for x in xs, y in ys]

heatmap(xs, ys, vals, xlabel="Re(z)", ylabel="Im(z)", title="|v(a,b,z,τ)|")


a, b = 0.3, 0.2
z = 0.2
N = 50
ys = range(0.1, 2, length=200)   # Im(τ)
vals = [abs(v(a,b,z,im*y,N)) for y in ys]

plot(ys, vals, xlabel="Im(τ)", ylabel="|v|")

z = 0.2
τ = 0.5im
N = 50
as = range(0, 1, length=100)
bs = range(0, 1, length=100)

vals = [abs(v(a,b,z,τ,N)) for a in as, b in bs]

heatmap(as, bs, vals, xlabel="a", ylabel="b", title="|v(a,b,z,τ)|")
