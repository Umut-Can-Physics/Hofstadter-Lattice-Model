"""
Get occupied particle site indices from many-body basis for a given basis index 
"""
ptl_site_ind(mb_basis, bi) = findall(x -> x≠0, mb_basis.occupations[bi])

"""
Calculate sum of y-coordinates squared for occupied sites in a given many-body basis index
"""
∑(Lattice, mb_basis, bi) = sum( Lattice.coordinates[:,2][ptl_site_ind(mb_basis,bi)] .^2 )

"""
Calculate sum of x-coordinates for occupied sites in a given many-body basis index
"""
∑x(Lattice, mb_basis, bi) = sum( Lattice.coordinates[:,1][ptl_site_ind(mb_basis,bi)])

"""
Calculate exponential factor for occupied sites in a given many-body basis index
"""
e(Lattice, mb_basis, bi, lb) = exp( - ∑(Lattice, mb_basis, bi) / (2*lb^2) ) 

"""   
Calculate exponential factor of the Composite Boson wavefunction for occupied sites in a given many-body basis index
"""
e_CB(Lattice, lb, lb_prime, bi) = exp( - ∑(Lattice, mb_basis, bi)/2 * ( 1/lb^2 - 1/lb_prime^2 ) )

"""
Calculate sum of z-coordinates for occupied sites in a given many-body basis index
"""
function Z_cm(Lattice, mb_basis, bi)
    return sum(Lattice.z_coords[ptl_site_ind(mb_basis,bi)])
end

"""
Calculate single-particle wave function in terms of Jacobi theta function
"""
function SpAnalyticWaveFunction(k, Nx, Ny, lat, UpperLimit, Nphi_prime, τ)
    return exp.(-pi * Nphi_prime/(Nx*Ny) .* (imag.(lat.z_coords)).^2).*v.(k/Nphi_prime, 0, Nphi_prime*lat.z_coords/Nx, τ*Nphi_prime, UpperLimit)
end

"""
Calculate trial center of mass wave function in terms of Jacobi theta function
"""
function CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, ar)
    Zcm = Z_cm(Lattice, mb_basis, bi)
    return v(ar, 0, 2*Zcm / Nx, 2*im*Ny/Nx, 20)*exp(-0*(4*pi*im*ar/Nx)*real.(Zcm))
end

"""
Calculate trial center of mass wave function in terms of Jacobi theta function
"""
function CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, ar)
    Zcm = Z_cm(Lattice, mb_basis, bi)
    return v(0, -ar/2, 2*(Zcm/Nx + ar/2), 2*im*Ny/Nx, 20)
end

"""
Calculate Composite Boson wave function for a given many-body basis index (components)
"""
function CB_Component(Lattice, mb_basis, bi, 𝜓ₛₚ, UpperLimit, lb, lb_prime, qq, pp)
    #- SP PART -#

    # find lattice index of these coordinates
    ptl_site_ind_vecs = ptl_site_ind(mb_basis, bi)
    # coordinates of occupied sites in the MB basis
    r1 = Lattice.coordinates[ptl_site_ind_vecs[1]]
    r2 = Lattice.coordinates[ptl_site_ind_vecs[2]]
    r = [r1, r2]
    
    # SP basis amplitudes corresponding to occupied MB basis
    α₁ = 𝜓ₛₚ[:,1][ptl_site_ind_vecs]
    α₂ = 𝜓ₛₚ[:,2][ptl_site_ind_vecs]

    𝜓₁₁ = α₁[1] * α₁[2] # 2 \phi_1(r_1) * \phi_1(r_2) 
    𝜓₂₂ = α₂[1] * α₂[2] # 2 \phi_2(r_1) * \phi_2(r_2) 
    𝜓₁₂ = α₁[1] * α₂[2] + α₂[1] * α₁[2] # \phi_1(r_1) * \phi_2(r_2) + \phi_2(r_1) * \phi_1(r_2)
    # 1/√N! * permanent(ϕ(r_1,r_2,..,r_N))

    𝜓ᵢⱼ = [𝜓₁₁, 𝜓₂₂, 𝜓₁₂] 

    # MB PART #

    ψ_rel = Relative(pn, ptl_site_ind(mb_basis,bi), lat, UpperLimit)
  
    WF = "CB"

    ψ_CM_list = [

    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 0/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 2/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 4/6),

    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 0/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 2/6)*exp(-2*pi*im/3)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 4/6)*exp(-4*pi*im/3),

    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 0/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 2/6)*exp(-4*pi*im/3)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 4/6)*exp(-8*pi*im/3), 
    
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 1/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 3/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 5/6),

    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 1/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 3/6)*exp(-2*pi*im/3)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 5/6)*exp(-4*pi*im/3),

    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 1/6)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 3/6)*exp(-4*pi*im/3)+
    CMAnsatz(Lattice, mb_basis, bi, Nx, Ny, 5/6)*exp(-8*pi*im/3), 

    #= CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, 0/6),
    CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, 1/6),
    CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, 2/6),
    CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, 3/6),
    CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, 4/6),
    CMAnsatz2(Lattice, mb_basis, bi, Nx, Ny, 5/6) =#

    ]

    ExpFun = e_CB(Lattice, lb, lb_prime, bi)

    ψ_mb = [ ψ_rel*ψ_CM_list[i]*ExpFun for i in eachindex(ψ_CM_list)]

    return kron(𝜓ᵢⱼ, ψ_mb) # \phi_1(r_1) * \phi_1(r_2) * \psi_L(r_1,r_2)
end

"""
Construct the ansatz wavefunction vector for given many-body basis and parameters
"""
function Ansatz(Lattice, mb_basis, UpperLimit, lb, lb_prime, WF, 𝜓ₛₚ, qq, number_of_ansatz, pp) 

    if WF == "Laughlin" 
        
        ψ_rel = [Relative(pn, ptl_site_ind(mb_basis,bi), lat, UpperLimit) for bi in 1:length(mb_basis)]
 
        ExpFun = [e(lat, mb_basis, bi, param.lb) for bi in 1:length(mb_basis)]

        CM0 = [CenterOfMass(mb_basis, bi, lat, param.Nphi, 0, qq, UpperLimit, WF) for bi in 1:length(mb_basis)]

        ψ0 = normalize(ψ_rel.*ExpFun.*CM0)
 
        CM1 = [CenterOfMass(mb_basis, bi, lat, param.Nphi, 1, qq, UpperLimit, WF) for bi in 1:length(mb_basis)]

        ψ1 = normalize(ψ_rel.*ExpFun.*CM1)

        Ansatz = [ψ0 ψ1] # two degenerate ground state at ν=1/2

    elseif WF == "CB"
        #@warn "CompositeBoson function is valid only for 2 particle and two sp ground state and degeneracy."

        Ansatz = zeros(ComplexF64, length(mb_basis), number_of_ansatz) # 3(sp) * 3(MB)

        for bi in eachindex(mb_basis.occupations)   
            Ansatz[bi,:] = CB_Component(Lattice, mb_basis, bi, 𝜓ₛₚ, UpperLimit, lb, lb_prime, qq, pp)
        end 
    end
    
    return Ansatz
end