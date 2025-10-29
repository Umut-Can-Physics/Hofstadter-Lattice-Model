"""
The definition of the Jacobi-Theta function
"""
function v(a, b, z, τ, UpperLimit)
    ThetaFun = [exp(im*pi*τ*(n+a).^2 + 2*pi*im*(n+a)*(z+b)) for n in -UpperLimit:UpperLimit] 
    return sum(ThetaFun)    
end

"""
Relative part of the generalized Laughlin wave function
"""
function Relative(pn, occupied_site_ind, Lattice, UpperLimit)
    a = b = 1/2
    τ = im*Lattice.Ny/Lattice.Nx
    Π = 1  
    # i<j
    for i in 1:pn
        for j in (i+1):pn 
            # zi - zj 
            dz = Lattice.z_coords[occupied_site_ind[i]] - Lattice.z_coords[occupied_site_ind[j]]
            Π *= v(a, b, dz/Lattice.Nx, τ, UpperLimit)^2
        end
    end 
    return Π
end 

"""
Center of mass part of the generalized Laughlin and CB wave functions
"""
function CenterOfMass(mb_basis, bi, Lattice, Nphi, d, qq, UpperLimit, WF)
    θ = 0
    if WF == "Laughlin"    
        a = d/2 + (Nphi-2)/4 # d=0 or d=1 refer to two degenerate ground state at ν=1/2
        b = -(Nphi-2)/2 
        z = 2 * sum(Lattice.z_coords[ptl_site_ind(mb_basis,bi)]) / Nx
        τ = 2*im*Ny/Nx
        θ = v(a, b, z, τ, UpperLimit) 
    elseif WF == "CB"
        a = d/qq + (Nphi-qq)/(2*qq)
        b = -(Nphi-qq)/2  
        z = qq * sum(Lattice.z_coords[ptl_site_ind(mb_basis,bi)]) / Nx
        τ = qq*im*Ny/Nx
        θ = v(a, b, z, τ, UpperLimit)  
    end
    return θ
end