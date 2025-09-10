# The definition of the Jacobi-Theta function
function v(a, b, z, τ, UpperLimit)
    ThetaFun = [exp(im*pi*τ*(n+a).^2 + 2*pi*im*(n+a)*(z+b)) for n in -UpperLimit:UpperLimit] 
    return sum(ThetaFun) 
end

# Sum of complex coordinates of all particles
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
function CenterOfMass(basis, bi, Nx, Ny, Nphi, d, α, UpperLimit, shift_amount, type, WF)
    θ = 0
    if WF == "Laughlin"    
        a = d/2 + (Nphi-2)/4 # d=0 or d=1 refer to two degenerate ground state at ν=1/2
        b = -(Nphi-2)/2 
        z = 2 * Z(basis, bi, type, shift_amount) / Nx
        τ = 2*im*Ny/Nx
        θ = v(a, b, z, τ, UpperLimit) 
    elseif WF == "CB"
        #qq = denominator(rationalize(pn/Nphi)) # Denominator of the filling factor ν
        qq = 2
        a = d/qq + (Nphi-qq)/(2*qq)
        b = -(Nphi-qq)/2 
        z = qq * Z(basis, bi, type, shift_amount) / Nx
        τ = qq*im*Ny/Nx
        θ = v(a, b, z, τ, UpperLimit) 
    end
    return θ
end