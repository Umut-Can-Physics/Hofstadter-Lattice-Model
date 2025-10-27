using Plots
using LaTeXStrings
includet("Main Scripts/Hofstadter.jl")
using .Hofstadter
include("Run.jl") # Build a model

UpperLimit = 10
WF ="Laughlin"
ansatz = Ansatz(lat, mb_basis, UpperLimit, param.lb, 0, WF, nothing)
ψ0 = ansatz[:,1]
ψ1 = ansatz[:,2]

Nev = Int(param.GroundStateDegeneracy) + 6
method = "Arpack" # "Lapack", "Arpack", "KrylovKit"
ϵ, ψ = SolveMatrix(HH, Nev, method)
scatter(real(ϵ))

# Note: Solve the model using Lapack, or use Arpack and optimize coefficient of eigenstates.

ψ0'*ψ1

scatter(abs.(ψ0), label=L"\psi_{d=0}", xlabel="Basis order", ylabel=L"|\psi|")
scatter!(abs.(ψ1))

scatter!(abs.(ψ[:,1]),label=L"\psi_{ED,1}", xlabel="Basis order", ylabel=L"|\psi|")
scatter!(abs.(ψ[:,2]),label=L"\psi_{ED,2}", xlabel="Basis order", ylabel=L"|\psi|")

heatmap(RealSpaceDensity(Nx, Ny, ψ0, OccBasis).+RealSpaceDensity(Nx, Ny, ψ1, OccBasis))

heatmap(RealSpaceDensity(Nx, Ny, ψ[:,1], OccBasis).+RealSpaceDensity(Nx, Ny, ψ[:,2], OccBasis))

Overlap(ψ0, ψ[:,1])
Overlap(ψ0, ψ[:,2])
Overlap(ψ1, ψ[:,1])
Overlap(ψ1, ψ[:,2])

W = OverlapMat(ψ0, ψ1, ψ[:,1], ψ[:,2])
HilbertSchmidtNorm(W)

overlap_values = CoeffOptimization(ψ0, ψ1, ψ[:,2])
maximum(overlap_values)