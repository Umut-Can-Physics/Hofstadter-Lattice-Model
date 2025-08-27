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

# Just convert matrix to vector
SiteCoordinates(Nx::Int, Ny::Int) = [Vector{Float64}(row) for row in eachrow(square_lattice(Nx,Ny)[2])]

function ParCoord(Nx::Int, Ny::Int, basis, i::Int, type::String)
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

e(Nx, Ny, lb, basis, bi, type, shift_amount) = exp( - sum( imag.( ComplexCoords(Nx, Ny, basis, bi, type, shift_amount) ).^2 ) / (2*lb^2) ) 

function GeneralizedLaughlin(basis, Nx, Ny, UpperLimit, type)
    
    ψ_rel = [Relative(pn, basis, bi, Nx, Ny, type, UpperLimit, shift_amount) for bi in 1:length(basis)]

    ExpFun = [e(Nx, Ny, lb, basis, i, type, shift_amount) for i in 1:length(basis)]

    d = 0
    ψ_CM0 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ0 = normalize(ψ_rel.*ψ_CM0.*ExpFun)
 
    d = 1
    ψ_CM1 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ1 = normalize(ψ_rel.*ψ_CM1.*ExpFun)

    return ψ0, ψ1
end

function CompositeBosonMBPart(basis, Nx, Ny, UpperLimit, type)
    ψ_rel = [Relative(pn, basis, bi, Nx, Ny, type, UpperLimit, shift_amount) for bi in 1:length(basis)]

    d = 0
    ψ_CM0 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ0 = normalize(ψ_rel.*ψ_CM0)
 
    d = 1
    ψ_CM1 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ1 = normalize(ψ_rel.*ψ_CM1)

    d = 2
    ψ_CM2 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ2 = normalize(ψ_rel.*ψ_CM2.*ExpFun)

    return ψ0, ψ1, ψ2
end