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

# Solve Hubbard-Hofstadter Model
#E, ψ = Solve(pn, Nx, Ny, alpha, periodicity, HardCore, U, Nev)