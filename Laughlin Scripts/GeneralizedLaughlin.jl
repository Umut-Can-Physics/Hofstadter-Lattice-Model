using Combinatorics

function BosonLabels(A::Vector{<:Real})
    result = Int[]
    for (idx, val) in enumerate(A)
        if val > 0
            num_repetitions = Int(val)
            for _ in 1:num_repetitions
                push!(result, idx)
            end
        end
    end
    return result
end

# Origion is put on the center of the lattice since the quasihole coordinate is 0+0i
# Note that Nx and Ny should be odd, because coordinates start from the center of the lattice.
# SiteCoordinates(Nx, Ny) = [ [x,y] for x in -(Nx-1)/2:(Nx-1)/2 for y in -(Ny-1)/2:(Ny-1)/2 ]
# Just convert matrix to vector
SiteCoordinates(Nx, Ny) = [Vector{Float64}(row) for row in eachrow(square_lattice(Nx,Ny)[2])]

function ParCoord(Nx, Ny, basis, i, type)
    SiteCoords = SiteCoordinates(Nx, Ny)  
    if type == "boson"
        result = SiteCoords[BosonLabels(basis[i])]
    elseif type=="fermion"
        result = SiteCoords[findall(x->x≠0, basis[i])]
    end
    return result
end

ComplexCoords(Nx, Ny, basis, i, type) = Complex.(getindex.(ParCoord(Nx, Ny, basis, i, type), 1),getindex.(ParCoord(Nx, Ny, basis, i, type), 2)).-(Nx/2-(1/2)+im*(Ny/2-1/2))

e(Nx, Ny, lb, basis, i, type) = exp(-sum(imag.(ComplexCoords(Nx, Ny, basis, i, type))).^2) /(2*lb^2)

# combination function automatically use the i<j condition
Comb(Nx, Ny, basis, i, type) = collect(combinations(ComplexCoords(Nx, Ny, basis, i, type), 2))

# z_i - z_j pairs such that i<j
ComplexCoordsDiff(Nx, Ny, basis, i, type) = getindex.(Comb(Nx, Ny, basis, i, type), 1) - getindex.(Comb(Nx, Ny, basis, i, type), 2)

function GeneralizedLaughlin(basis, Nx, Ny, UpperLimit, type)
    
    ψ_rel = [Relative(basis, i, Nx, Ny, UpperLimit, type) for i in 1:length(basis)]

    ExpFun = [e(Nx, Ny, lb, basis, i, type) for i in 1:length(basis)]

    l = 0
    ψ_CM0 = [CenterOfMass(basis, i, Nx, Ny, l, alpha, UpperLimit, type) for i in 1:length(basis)]
    ψ0 = normalize(ψ_rel.*only.(ψ_CM0).*ExpFun)
 
    l = 1
    ψ_CM1 = [CenterOfMass(basis, i, Nx, Ny, l, alpha, UpperLimit, type).*e(Nx, Ny, lb, basis, i, type) for i in 1:length(basis)]
    ψ1 = normalize(ψ_rel.*only.(ψ_CM1).*ExpFun)
    
    return ψ0, ψ1
end