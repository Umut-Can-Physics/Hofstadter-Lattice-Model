using Revise
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
periodicity = true
HardCore = true
Nev = 10
gauge = "Landau" # "Landau", "Symmetric"
imp_str = 0.1
perturbation = true

# Solve Hubbard-Hofstadter Model
method = "Lapack" # "Lapack", "Arpack", "KrylovKit"
E_Lapack, ψ_Lapack = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method)

method = "Arpack" 
E_Arpack, ψ_Arpack = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method)

method = "KrylovKit"
E_KrylovKit, ψ_KrylovKit = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method)

scatter(E_Lapack[1:Nev], label="Lapack")
scatter!(real.(E_Arpack), label="Arpack", markershape=:star)
scatter!(E_KrylovKit, label="KrylovKit", markershape=:diamond)

ψ_KrylovKit_Matrix = hcat([ψ_KrylovKit[i].data for i in 1:size(ψ_KrylovKit, 1)]...)

scatter(abs.(ψ_Lapack[:,1]), label="Lapack", xlabel="Basis order", ylabel="|ψ|", title="Two Fold Degeneracy d=1, perturbation=$perturbation")
scatter!(abs.(ψ_Arpack[:,1]), label="Arpack", markershape=:star)
scatter!(abs.(ψ_KrylovKit_Matrix[:,1]), label="KrylovKit", markershape=:diamond)

scatter(abs.(ψ_Lapack[:,2]), label="Lapack", xlabel="Basis order", ylabel="|ψ|", title="Two Fold Degeneracy d=2, perturbation=$perturbation")
scatter!(abs.(ψ_Arpack[:,2]), label="Arpack", markershape=:star)
scatter!(abs.(ψ_KrylovKit_Matrix[:,2]), label="KrylovKit", markershape=:diamond)