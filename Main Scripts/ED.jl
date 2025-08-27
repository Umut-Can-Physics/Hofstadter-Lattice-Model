using LinearAlgebra
using Arpack

"""
Sort eigenstates by real energies
"""
function SortStates(E, ψ)
    try
        i = sortperm(E, by=real)
        E = E[i]
        ψ = ψ[:, i] 
    catch e
        # If a BoundsError occurs, catch it
        if isa(e, BoundsError)
            i = sortperm(E, by=real)
            E = E[i]
            ψ = ψ[i]
            return E, ψ
        else
            # Re-throw any other unexpected errors
            throw(e)
        end
    end
    return E, ψ
end

function Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
    if problem_type == "MB"
        H = HubbardHofstadter(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, perturbation, imp_str)
        if method == "Lapack"
            H = Matrix(H.data)
            E, ψ = eigen(H); 
            E, ψ = SortStates(E, ψ)
        elseif method == "Arpack"
            H = Matrix(H.data)
            E, ψ = eigs(H, nev=Nev, which=:SR)
            E, ψ = SortStates(E, ψ)
        elseif method == "KrylovKit"
            E, ψ = eigenstates(dense(H), Nev)
            E, ψ = SortStates(E, ψ)
        else
            error("Invalid method. Use 'Lapack', 'Arpack', or 'KrylovKit'.")
        end
    elseif problem_type == "SP"
        H = SingleParticleModel(Nx, Ny, alpha, periodicity, gauge)
        E, ψ = eigs(H, nev=Nev, which=:SR)
        E, ψ = SortStates(E, ψ)
    end
    return E, ψ
end