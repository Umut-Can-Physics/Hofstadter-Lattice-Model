using Revise
using QuantumOptics
using Plots
using LaTeXStrings
includet("Laughlin Scripts/GeneralizedLaughlin.jl")
includet("Laughlin Scripts/JacobiThetaFunction.jl")
includet("Main Scripts/Lattice.jl")
includet("Main Scripts/Model.jl")
includet("Main Scripts/MBBasis.jl")
includet("Main Scripts/Operators.jl")
includet("Main Scripts/ED.jl")
includet("Main Scripts/Utilities.jl")

pn = 2
Nx = 4
Ny = 4
alpha = 1/4
U = 1
lb = 1/sqrt(2*pi*alpha)
periodicity = true
HardCore = true
Nev = 10
shift_amount = 0

imp_str = 0
perturbation = false
E, ψ = Solve(pn, Nx, Ny, alpha, periodicity, HardCore, U, Nev, perturbation, imp_str)

N = Nx*Ny
spbasis = NLevelBasis(N)
basis = fermionstates(spbasis, pn)
type = "fermion"
rel = []
cm = []
UpperLimit = 10

ψ0, ψ1 = GeneralizedLaughlin(basis, Nx, Ny, UpperLimit, type)

ψ0'*ψ1 #check

scatter(abs.(ψ0), label=L"\psi_{d=0}", xlabel="Basis order", ylabel=L"|\psi|")
scatter!(abs.(ψ1))

scatter!(abs.(ψ[:,1]),label=L"\psi_{ED,1}", xlabel="Basis order", ylabel=L"|\psi|")
scatter!(abs.(ψ[:,2]),label=L"\psi_{ED,2}", xlabel="Basis order", ylabel=L"|\psi|")

mb = MBBasis(pn, Nx, Ny, HardCore)

heatmap(RealSpaceDensity(Nx, Ny, ψ0, mb).+RealSpaceDensity(Nx, Ny, ψ1, mb))

heatmap(RealSpaceDensity(Nx, Ny, ψ[:,1], mb).+RealSpaceDensity(Nx, Ny, ψ[:,2], mb))

Overlap(ψ0, ψ[:,1])
Overlap(ψ0, ψ[:,2])
Overlap(ψ1, ψ[:,1])
Overlap(ψ1, ψ[:,2])

W = OverlapMat(ψ0, ψ1, ψ[:,1], ψ[:,2])
HilbertSchmidtNorm(W)

overlap_values = CoeffOptimization(ψ0, ψ1, ψ[:,2])
maximum(overlap_values)
plot(overlap_values)