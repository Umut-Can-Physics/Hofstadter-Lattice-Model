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
    return abs.(ψ1'*ψ2)^2
end

"Overlap matrix for Laughlin ν=1/2"
function OverlapMat(ψ0, ψ1, ψED_1, ψED_2)
    # Generalized Laughlin Wave Functions for d=0 and d=1
    AnalyticWaveFunction = [ψ0, ψ1]
    # ED States
    EDStates = [ψED_1, ψED_2]
    OvMat = zeros(2,2) 
    for i in 1:2
        for j in 1:2
            OvMat[i,j] = Overlap(AnalyticWaveFunction[i], EDStates[j])
        end
    end
    return OvMat
end

function HilbertSchmidtNorm(W)
    dim = size(W, 1) # Dimension of square matrix
    Σ = 0
    for i in 1:dim
        for j in 1:dim
            Σ += abs((W'*W)[i,j])^2
        end
    end
    return sqrt((1/dim)*Σ)
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