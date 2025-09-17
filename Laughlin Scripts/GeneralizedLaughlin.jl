using Combinatorics

ptl_site_ind(mb_basis,bi) = findall(x->x≠0, mb_basis.occupations[bi])

∑(Lattice, mb_basis, bi) = sum( Lattice.coordinates[:,2][ptl_site_ind(mb_basis,bi)] .^2 )

e(Lattice, mb_basis, bi, lb) = exp( - ∑(Lattice, mb_basis, bi) / (2*lb^2) ) 

e_CB(Lattice, lb, lb_prime, bi) = exp( - ∑(Lattice, mb_basis, bi)/2 * ( 1/lb^2 - 1/lb_prime^2 ) )

function CB_Component(Lattice, mb_basis, bi, 𝜓ₛₚ, UpperLimit, lb, lb_prime)
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
    𝜓₁₁ = 1 * α₁[1] * α₁[2] # 2 \phi_1(r_1) * \phi_1(r_2) 
    𝜓₂₂ = 1 * α₂[1] * α₂[2] # 2 \phi_2(r_1) * \phi_2(r_2) 
    𝜓₁₂ = (α₁[1] * α₂[2] + α₂[1] * α₁[2])/sqrt(2) # \phi_1(r_1) * \phi_2(r_2) + \phi_2(r_1) * \phi_1(r_2)
    # 1/√N! * permanent(ϕ(r_1,r_2,..,r_N))

    𝜓ᵢⱼ = [𝜓₁₁, 𝜓₂₂, 𝜓₁₂] 

    # MB PART #

    ψ_rel = Relative(pn, ptl_site_ind(mb_basis,bi), lat, UpperLimit)
  
    WF = "CB"
    ψ_CM0 = CenterOfMass(mb_basis, bi, lat, param.Nphi, 0, UpperLimit, WF)
    ψ_CM1 = CenterOfMass(mb_basis, bi, lat, param.Nphi, 1, UpperLimit, WF)
    #ψ_CM2 = CenterOfMass(mb_basis, bi, lat, param.Nphi, 2, UpperLimit, WF)
    Zcm = sum(Lattice.z_coords[ptl_site_ind(mb_basis,bi)])
    ψ_CM2 = v(1/2, -1/2, Zcm/Nx, im*Ny/Nx, UpperLimit)   

    ExpFun = e_CB(Lattice, lb, lb_prime, bi)

    𝜓₁ = ψ_rel*ψ_CM0*ExpFun # \psi_L(r_1,r_2) for d=0
    𝜓₂ = ψ_rel*ψ_CM1*ExpFun # \psi_L(r_1,r_2) for d=1
    𝜓₃ = ψ_rel*ψ_CM2*ExpFun # \psi_L(r_1,r_2) for d=2

    𝜓ₐ = [𝜓₁, 𝜓₂, 𝜓₃] 

    return kron(𝜓ᵢⱼ, 𝜓ₐ) # \phi_1(r_1) * \phi_1(r_2) * \psi_L(r_1,r_2)
end

function Ansatz(Lattice, mb_basis, UpperLimit, lb, lb_prime, WF, 𝜓ₛₚ)

    if WF == "Laughlin" 
        
        ψ_rel = [Relative(pn, ptl_site_ind(mb_basis,bi), lat, UpperLimit) for bi in 1:length(mb_basis)]
 
        ExpFun = [e(lat, mb_basis, bi, param.lb) for bi in 1:length(mb_basis)]

        CM0 = [CenterOfMass(mb_basis, bi, lat, param.Nphi, 0, UpperLimit, WF) for bi in 1:length(mb_basis)]

        ψ0 = normalize(ψ_rel.*ExpFun.*CM0)
 
        CM1 = [CenterOfMass(mb_basis, bi, lat, param.Nphi, 1, UpperLimit, WF) for bi in 1:length(mb_basis)]

        ψ1 = normalize(ψ_rel.*ExpFun.*CM1)

        Ansatz = [ψ0 ψ1] # two degenerate ground state at ν=1/2

    elseif WF == "CB"
        @warn "CompositeBoson function is valid only for 2 particle and two sp ground state and degeneracy."

        Ansatz = zeros(ComplexF64, length(mb_basis), 9) # 3(sp) * 3(MB)

        for bi in eachindex(mb_basis.occupations)
            Ansatz[bi,:] = CB_Component(Lattice, mb_basis, bi, 𝜓ₛₚ, UpperLimit, lb, lb_prime)
        end 
    end
    return Ansatz
end