using LinearAlgebra

function Solve(pn, Nx, Ny, alpha, periodicity, HardCore, U, Nev)
    H = HofstadterHubbard(pn, Nx, Ny, alpha, periodicity, HardCore, U)
    E, ψ = eigenstates(dense(H), Nev); 
    return E, ψ
end