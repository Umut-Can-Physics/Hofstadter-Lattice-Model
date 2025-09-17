using LinearAlgebra
using Arpack

"""
Sort eigenstates by real energies
"""
function SortStates(E, ψ)
    i = sortperm(E, by=real)
    E = E[i]
    if length(size(ψ)) == 2 # If data is matrix object
        ψ = ψ[:, i] 
    elseif length(size(ψ)) == 1 # If data is QoJulia object
        ψ = ψ[i]
    end

    return E, ψ
end

function SolveMatrix(H, Nev, method)
    if method == "Lapack"
        println("Hi!")
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
    return E, ψ
end