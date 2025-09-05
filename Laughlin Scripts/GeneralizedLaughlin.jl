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

∑(Nx, Ny, basis, bi, type, shift_amount) = sum( imag.( ComplexCoords(Nx, Ny, basis, bi, type, shift_amount) ).^2 )
e_CB(Nx, Ny, lb, lb_prime, basis, bi, type, shift_amount) = exp( - ∑(Nx, Ny, basis, bi, type, shift_amount)/2 * ( 1/lb^2 - 1/lb_prime^2 ) )

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

    ExpFun = [e_CB(Nx, Ny, lb, lb_prime, basis, bi, type, shift_amount) for bi in 1:length(basis)]
 
    d = 0
    ψ_CM0 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ0 = ψ_rel.*ψ_CM0.*ExpFun
 
    d = 1
    ψ_CM1 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ1 = ψ_rel.*ψ_CM1.*ExpFun

    d = 2
    ψ_CM2 = [CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type) for bi in 1:length(basis)]
    ψ2 = ψ_rel.*ψ_CM2.*ExpFun

    return ψ0, ψ1, ψ2
end

# 𝜓ₛₚ two degenerate ground state
function CompositeBoson(basis, Nx, Ny, lb, lb_prime, UpperLimit, type, 𝜓ₛₚ, z_nu_row0, z_nu_row1, z_nu_row2)

    #@warn "CompositeBoson function is valid only for 2 particle and two sp ground state and degeneracy."
    
    SiteCoords = SiteCoordinates(Nx, Ny)

    ψ_CB = zeros(ComplexF64, length(basis), 9) # 3(sp) * 3(MB)

    for bi in eachindex(basis)

        #- SP PART -#

        # coordinates of occupied sites in the MB basis
        r1 = [real(ComplexCoords(Nx, Ny, basis, bi, type, shift_amount)[1]), imag(ComplexCoords(Nx, Ny, basis, bi, type, shift_amount)[1])]
        r2 = [real(ComplexCoords(Nx, Ny, basis, bi, type, shift_amount)[2]), imag(ComplexCoords(Nx, Ny, basis, bi, type, shift_amount)[2])]
        r = [r1, r2]
        # find lattice index of these coordinates
        R = indexin(r, SiteCoords) 
        
        # SP basis amplitudes corresponding to occupied MB basis
        α₁ = 𝜓ₛₚ[:,1][R]
        α₂ = 𝜓ₛₚ[:,2][R]
        𝜓₁₁ = α₁[1] * α₁[2] # \phi_1(r_1) * \phi_1(r_2) 
        𝜓₂₂ = α₂[1] * α₂[2] # \phi_2(r_1) * \phi_2(r_2) 
        𝜓₁₂ = α₁[1] * α₂[2] # \phi_1(r_1) * \phi_2(r_2) 

        𝜓ᵢⱼ = [𝜓₁₁, 𝜓₂₂, 𝜓₁₂]

        # MB PART #
 
        ψ_rel = Relative(pn, basis, bi, Nx, Ny, type, UpperLimit, shift_amount)

        
        #= ψ_CM0 = CenterOfMass(basis, bi, Nx, Ny, 0, alpha, UpperLimit, shift_amount, type)
        ψ_CM1 = CenterOfMass(basis, bi, Nx, Ny, 1, alpha, UpperLimit, shift_amount, type)
        ψ_CM2 = CenterOfMass(basis, bi, Nx, Ny, 2, alpha, UpperLimit, shift_amount, type) =#
        

        #= ψ_CM0 = CenterOfMass_CB(basis, bi, Nx, Ny, pn, Nphi, 0, alpha, UpperLimit, shift_amount, type)
        ψ_CM1 = CenterOfMass_CB(basis, bi, Nx, Ny, pn, Nphi, 1, alpha, UpperLimit, shift_amount, type)
        ψ_CM2 = CenterOfMass_CB(basis, bi, Nx, Ny, pn, Nphi, 2, alpha, UpperLimit, shift_amount, type) =#

        ψ_CM0 = CM_New_2(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row0)
        ψ_CM1 = CM_New_2(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row1)
        ψ_CM2 = CM_New_2(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row2)

        #= ψ_CM0 = CM_New(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row0)
        ψ_CM1 = CM_New(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row1)
        ψ_CM2 = CM_New(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row2) =#

        ExpFun = e_CB(Nx, Ny, lb, lb_prime, basis, bi, type, shift_amount)

        𝜓₁ = ψ_rel*ψ_CM0*ExpFun # \psi_L(r_1,r_2) for d=0
        𝜓₂ = ψ_rel*ψ_CM1*ExpFun # \psi_L(r_1,r_2) for d=1
        𝜓₃ = ψ_rel*ψ_CM2*ExpFun # \psi_L(r_1,r_2) for d=2

        𝜓ₐ = [𝜓₁, 𝜓₂, 𝜓₃]

        ψ_CB[bi,:] = kron(𝜓ᵢⱼ, 𝜓ₐ) # \phi_1(r_1) * \phi_1(r_2) * \psi_L(r_1,r_2)

    end

    return ψ_CB
end