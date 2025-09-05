# The definition of the Jacobi-Theta function
# 
function v(a, b, z, τ, UpperLimit)
    ThetaFun = [exp(im*pi*τ*(n+a).^2 + 2*pi*im*(n+a)*(z+b)) for n in -UpperLimit:UpperLimit]
    return sum(ThetaFun) 
end
 
function Z(basis, i, type, shift_amount)
    return sum(ComplexCoords(Nx, Ny, basis, i, type, shift_amount))
end

# Relative part of the generalized Laughlin wave function
function Relative(pn, basis, bi, Nx, Ny, type, UpperLimit, shift_amount)
    a = b = 1/2
    τ = im*Ny/Nx 
    Π = 1 
    # i<j
    for i in 1:pn
        for j in (i+1):pn 
            # zi - zj 
            difference = ComplexCoords(Nx, Ny, basis, bi, type, shift_amount)[i] - ComplexCoords(Nx, Ny, basis, bi, type, shift_amount)[j]
            z = difference / Nx
            Π *= v(a, b, z, τ, UpperLimit)^2
        end
    end
    return Π
end 

# Center of mass part of the Laughlin wave function
# d=0 or d=1 refer to two degenerate ground state at \nu=1/2
function CenterOfMass(basis, bi, Nx, Ny, d, alpha, UpperLimit, shift_amount, type)
    N = Nx*Ny
    Nϕ = N*alpha
    a = d/2 + (Nϕ-2)/4
    b = -(Nϕ-2)/2 
    z = 2 * Z(basis, bi, type, shift_amount) / Nx
    τ = 2*im*Ny/Nx
    return v(a, b, z, τ, UpperLimit) 
end

function CenterOfMass_CB(basis, bi, Nx, Ny, pn, Nphi, d, alpha, UpperLimit, shift_amount, type)
    qq = denominator(rationalize(pn/Nphi)) 
    N = Nx*Ny 
    Nϕ = N*alpha
    a = d/qq + (Nϕ-qq)/(2*qq)
    b = -(Nϕ-qq)/2 
    z = qq * Z(basis, bi, type, shift_amount) / Nx
    τ = qq*im*Ny/Nx
    return v(a, b, z, τ, UpperLimit) 
end

function CM_New(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu_row)
    a = b = 1/2
    τ = im*Ny/Nx 
    Π = 1
    for i in 1:3
        z = (Z(basis, bi, type, shift_amount)-z_nu_row[i]) / Nx
        Π *= v(a, b, z, τ, UpperLimit)
    end
    return Π
end 

function CM_New_2(basis, bi, Nx, Ny, type, UpperLimit, shift_amount, z_nu)
    a = b = 1/2
    τ = im*Ny/Nx 
    z = ( Z(basis, bi, type, shift_amount) - z_nu ) / Nx
    return v(a, b, z, τ, UpperLimit)^2
end 

function CM_New_3(Nx, Ny, UpperLimit, z_coord, z_nu)
    a = b = 1/2
    τ = im*Ny/Nx 
    z = ( z_coord - z_nu ) / Nx
    return v(a, b, z, τ, UpperLimit)^2
end 