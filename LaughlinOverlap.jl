using Plots
using LaTeXStrings
using Serialization
includet("Main Scripts/Hofstadter.jl")
using .Hofstadter

include("Run.jl") # Build a sp model

if pn/param.Nphi == 1/2
    println("Running Laughlin Overlap Calculation at ν=1/2")    
else
    error("This Laughlin Overlap script is currently set up only for ν=1/2. Please adjust parameters in Run.jl accordingly.")
end

UpperLimit = 10
WF ="Laughlin"
ansatz = Ansatz(lat, mb_basis, UpperLimit, param.lb, nothing, WF, nothing, nothing, nothing, nothing)
ψ0 = ansatz[:,1]
ψ1 = ansatz[:,2]

Nev = Int(param.GroundStateDegeneracy) + 6
method = "Lapack" # "Lapack", "Arpack", "KrylovKit"
ϵ, ψ = SolveMatrix(HH, Nev, method)
scatter(real(ϵ))

# Note: Solve the model using Lapack, or use Arpack and optimize coefficient of eigenstates.

ψ0'*ψ1

scatter(abs.(ψ0), label=L"\psi_{d=0}", xlabel="Basis order", ylabel=L"|\psi|")
scatter!(abs.(ψ1))

scatter!(abs.(ψ[:,1]),label=L"\psi_{ED,1}", xlabel="Basis order", ylabel=L"|\psi|")
scatter!(abs.(ψ[:,2]),label=L"\psi_{ED,2}", xlabel="Basis order", ylabel=L"|\psi|")

heatmap(RealSpaceDensity(Nx, Ny, ψ0, mb_basis.occupations).+RealSpaceDensity(Nx, Ny, ψ1, mb_basis.occupations))

heatmap(RealSpaceDensity(Nx, Ny, ψ[:,1], OccBasis).+RealSpaceDensity(Nx, Ny, ψ[:,2], OccBasis))

Overlap(ψ0, ψ[:,1])
Overlap(ψ0, ψ[:,2])
Overlap(ψ1, ψ[:,1])
Overlap(ψ1, ψ[:,2])

W = OverlapMat(ψ0, ψ1, ψ[:,1], ψ[:,2])
HilbertSchmidtNorm(W)

overlap_values = CoeffOptimization(ψ0, ψ1, ψ[:,1])
maximum(overlap_values)