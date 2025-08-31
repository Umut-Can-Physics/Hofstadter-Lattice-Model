# The definition of the Jacobi-Theta function
function v(a, b, z, τ, UpperLimit)
    ThetaFun = [exp(im*pi*τ*(n+a).^2 + 2*pi*im*(n+a)*(z+b)) for n in -UpperLimit:UpperLimit]
    return sum(ThetaFun)
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

function Z(basis, i, type, shift_amount)
    return sum(ComplexCoords(Nx, Ny, basis, i, type, shift_amount))
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