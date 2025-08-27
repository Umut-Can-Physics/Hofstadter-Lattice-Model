include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots

pn = 2
Nx = 10
Ny = 10
alpha = -3/50
U = 1.0
periodicity = true
HardCore = true
Nphi = abs(Nx * Ny * alpha)
Nd = Int(Nphi-2*pn)
GroundStateDegeneracy = factorial(Nd + pn - 1) / (factorial(pn - 1)*factorial(Nd)) * (Nphi / pn)
Nev = Int(GroundStateDegeneracy) + 1
gauge = "Landau" # "Landau", "Symmetric"
imp_str = 0.1
perturbation = false

type = "fermion" # HardCore
UpperLimit = 10
shift_amount = 0
lb = 1/sqrt(2*pi*abs(alpha))
OccBasis = MBBasis(pn, Nx, Ny, HardCore).occupations
ψ0, ψ1, ψ2  = CompositeBosonMBPart(OccBasis, Nx, Ny, UpperLimit, type)

problem_type = "MB" # 'SP' or 'MB'
method = "Arpack" # "Lapack", "Arpack", "KrylovKit"
E_mb, ψ_mb = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_Lapack))

problem_type = "SP"
alpha = -1/50 
sp_basis = fermionstates(NLevelBasis(Nx*Ny), 1) # I suppose that this is the single particle basis consist of the single particle matrix
E_sp, ψ_sp = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_sp))

ψ_1 = ψ_sp[:,1]
ψ_2 = ψ_sp[:,2]

ψ_CB = ψ_1[z1] .* ψ_2[z2] .* ψ0 .* ψ1 .* ψ2