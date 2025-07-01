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

function ComplexCoords(Nx, Ny, basis, i, type, shift_amount)
    base_coords = Complex.(getindex.(ParCoord(Nx, Ny, basis, i, type), 1), getindex.(ParCoord(Nx, Ny, basis, i, type), 2))
    return base_coords .- shift_amount
end

e(Nx, Ny, lb, basis, i, type) = exp(-sum(imag.(ComplexCoords(Nx, Ny, basis, i, type, shift_amount))).^2) /(2*lb^2)

ComplexCoordsDiff(Nx, Ny, basis, i, type) = ComplexCoords(Nx, Ny, basis, i, type, shift_amount)

function GeneralizedLaughlin(basis, Nx, Ny, UpperLimit, type)
    
    ψ_rel = [Relative(basis, bi, Nx, Ny, type, shift_amount) for bi in 1:length(basis)]

    ExpFun = [e(Nx, Ny, lb, basis, i, type) for i in 1:length(basis)]

    l = 0
    ψ_CM0 = [CenterOfMass(basis, bi, Nx, Ny, l, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ0 = normalize(ψ_rel.*ψ_CM0.*ExpFun)

    l = 1
    ψ_CM1 = [CenterOfMass(basis, i, Nx, Ny, l, alpha, UpperLimit, shift_amount, type).*e(Nx, Ny, lb, basis, i, type) for i in 1:length(basis)]
    ψ1 = normalize(ψ_rel.*ψ_CM1.*ExpFun)
    
    return ψ0, ψ1
end