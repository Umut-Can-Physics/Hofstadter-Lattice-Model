# The definition of the Jacobi-Theta function
function v(a, b, z, τ, UpperLimit)
    Sum = zeros(length(z))
    for n in -UpperLimit:UpperLimit
        Sum = Sum .+ exp.(im*pi*τ*(n+a)^2 .+ 2*pi*im*(n+a).*(z.+b))
    end
    return Sum
end

# Relative part of the generalized Laughlin wave function
function Relative(basis, i, Nx, Ny, UpperLimit, type)
    a = b = 1/2
    z = ComplexCoordsDiff(Nx, Ny, basis, i, type)./Nx 
    τ = im*Ny/Nx
    return  prod( v(a, b, z, τ, UpperLimit).^2 ) 
end

function Z(basis, i, type)
    return sum( ComplexCoords(Nx, Ny, basis, i, type) )
end

# Center of mass part of the Laughlin wave function
# l=0 or l=1 refer to two degenerate ground state 
function CenterOfMass(basis, i, Nx, Ny, l, alpha, UpperLimit, type)
    N = Nx*Ny
    Nϕ = N*alpha
    a = l/2 + (Nϕ-2)/4
    b = -(Nϕ-2)/2
    z = 2 * Z(basis, i, type)/Nx
    τ = 2*im*Ny/Nx
    return v(a, b, z, τ, UpperLimit)
end