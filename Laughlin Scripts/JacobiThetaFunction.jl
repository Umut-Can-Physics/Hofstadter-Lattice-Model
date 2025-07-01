# The definition of the Jacobi-Theta function
function v(a, b, z, τ, UpperLimit)
    Sum = 0
    for n in -UpperLimit:UpperLimit
        Sum = Sum + exp(im*pi*τ*(n+a)^2 + 2*pi*im*(n+a)*(z+b))
    end
    return Sum
end

# Relative part of the generalized Laughlin wave function
function Relative(basis, bi, Nx, Ny, type, shift_amount)
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
# l=0 or l=1 refer to two degenerate ground state 
function CenterOfMass(basis, i, Nx, Ny, l, alpha, UpperLimit, shift_amount, type)
    N = Nx*Ny
    Nϕ = N*alpha
    a = l/2 + (Nϕ-2)/4
    b = -(Nϕ-2)/2
    z = 2 * Z(basis, i, type, shift_amount)/Nx
    τ = 2*im*Ny/Nx
    return v(a, b, z, τ, UpperLimit)
end