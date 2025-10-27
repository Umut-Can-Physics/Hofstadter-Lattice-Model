using LinearAlgebra

"""
Get matrix that have orthonormalized column vectors
"""
function QRDecomp(M::Matrix)
    F = qr(M)
    # The Q factor contains the orthonormalized columns
    Q = F.Q
    return Matrix(Q)
end

function Overlap(ψ1, ψ2)
    return ψ1'*ψ2
end

"Overlap matrix for Laughlin ν=1/2"
function OverlapMat(ψ0, ψ1, ψED_1, ψED_2)
    # Generalized Laughlin Wave Functions for d=0 and d=1
    AnalyticWaveFunction = [ψ0, ψ1]
    # ED States
    EDStates = [ψED_1, ψED_2]
    OvMat = zeros(Complex, 2,2) 
    for i in 1:2
        for j in 1:2
            OvMat[i,j] = Overlap(AnalyticWaveFunction[i], EDStates[j])
        end
    end
    return OvMat
end

function HilbertSchmidtNorm(W)
    dim = size(W, 1) # Dimension of square matrix
    WW = sum(W'*W)
    return (1/dim)*WW
end

function RealSpaceDensity(Nx, Ny, ψ,mb)
    N = Nx*Ny
    Density = zeros(N)
    for n in 1:N
        Density[n] = ψ'*number(mb, n).data*ψ
    end
    return reshape(Density, Nx, Ny)
end

"""
Optimize the coefficients to maximize the Laughlin overlap for two ground states degeneracies.
"""
function CoeffOptimization(ψ_1, ψ_2, Ψ_ED)
    overlap_values = []
    for r in 0:0.01:1
        for θ in 0:0.01:2*pi
            b = r*exp(im*θ)
            if abs(b) <= 1 # since a is real
                a = sqrt(1-abs(b)^2) 
                Ψ = a*ψ_1 + b*ψ_2
                Overlap = abs(Ψ'*Ψ_ED)^2
                push!(overlap_values, Overlap)
            end
        end
    end
    return overlap_values
end

function random_complex_triples(Nx::Int, Ny::Int, num_samples::Int)
    triples = []

    while length(triples) < num_samples
        # generate 3 random complex numbers
        z = Complex{Int}[]
        while length(z) < 3
            x, y = rand(0:Nx-1), rand(0:Ny-1)
            c = complex(x, y)
            if !(c in z) # ensure uniqueness
                push!(z, c)
            end
        end
        push!(triples, tuple(z...))
    end

    return triples
end

function gram_schmidt_qr(vectors::Vector{Vector{T}}) where T<:Number
    # Convert array-of-vectors into a matrix (columns are the vectors)
    M = hcat(vectors...)  
    
    # QR factorization
    F = qr(M)
    
    # Extract Q (orthonormal basis as a matrix)
    Q = Matrix(F.Q)
    
    # Return as array-of-vectors
    return [Q[:,i] for i in 1:size(Q,2)]
end