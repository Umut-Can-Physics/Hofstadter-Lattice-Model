using LinearAlgebra
using Arpack

function Solve(pn, Nx, Ny, alpha, periodicity, HardCore, U, Nev, perturbation, imp_str)
    H = HofstadterHubbard(pn, Nx, Ny, alpha, periodicity, HardCore, U, perturbation, imp_str)
    H = Matrix(H.data)
    E, ψ = eigen(H); 
    i = sortperm(E, by=real)
    E = E[i]
    ψ = ψ[:, i]
    return E, ψ
end