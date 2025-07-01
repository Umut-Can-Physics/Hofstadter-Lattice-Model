using Revise
using QuantumOptics
using Plots
includet("Laughlin Scripts/GeneralizedLaughlin.jl")
includet("Laughlin Scripts/JacobiThetaFunction.jl")
includet("Main Scripts/Lattice.jl")
includet("Main Scripts/Model.jl")
includet("Main Scripts/MBBasis.jl")
includet("Main Scripts/Operators.jl")
includet("Main Scripts/ED.jl")

pn = 2
Nx = 4
Ny = 4
alpha = 1/4
U = 1
lb = 1/sqrt(2*pi*alpha)
periodicity = true
HardCore = true
Nev = 10

E, ψ = Solve(pn, Nx, Ny, alpha, periodicity, HardCore, U, Nev)

function Overlap(ψ1, ψ2)
    return abs.(ψ1'*ψ2)^2
end

N = Nx*Ny
spbasis = NLevelBasis(N)
basis = fermionstates(spbasis, pn)
type = "fermion"
l = 0
rev = []
cm = []
UpperLimit = 100
for i in 1:length(basis)
    push!(rev, Relative(basis, i, Nx, Ny, UpperLimit, type))
    push!(cm, CenterOfMass(basis, i, Nx, Ny, l, alpha, UpperLimit, type))
end
scatter(only.(real(rev)))

ψ0, ψ1 = GeneralizedLaughlin(basis, Nx, Ny, UpperLimit, type)

mb = MBBasis(pn, Nx, Ny, HardCore)

function RealSpaceDensity(Nx, Ny, ψ,mb)
    N = Nx*Ny
    Density = zeros(N)
    for n in 1:N
        Density[n] = ψ'*number(mb, n).data*ψ
    end
    return reshape(Density, Nx, Ny)
end

heatmap(RealSpaceDensity(Nx, Ny, ψ0, mb))

heatmap(RealSpaceDensity(Nx, Ny, ψ1, mb))

heatmap(RealSpaceDensity(Nx, Ny, ψ[1].data, mb).+RealSpaceDensity(Nx, Ny, ψ[2].data, mb))

Overlap(ψ0, ψ[1].data)
Overlap(ψ0, ψ[2].data)
Overlap(ψ1, ψ[1].data)
Overlap(ψ1, ψ[2].data)