function SingleParticleModel(Lattice, α::Float64, gauge::String)
    
    if gauge != "Landau" && gauge != "Symmetric"
        error("Invalid gauge choice. Use 'Landau' or 'Symmetric'.")
    end

    @warn("Symmetric gauge is not implemented yet, use Landau gauge instead.")
 
    Nx = Lattice.Nx
    Ny = Lattice.Ny
    N = Lattice.N
    neig = Lattice.neig
    coordinates = Lattice.coordinates
     
    t = -1
    H = zeros(Complex{Float64},N,N)
    
    for m in 1:N
        for n in neig[m]

            # hopping from edge to edge
            if abs(coordinates[m,1]-coordinates[n,1])==Nx-1
                if coordinates[m,1] > coordinates[n,1]
                    H[m,n] = t*exp(-1im*2*pi*α*coordinates[m,2])
                elseif coordinates[m,1] < coordinates[n,1]
                    H[m,n] = t*exp(1im*2*pi*α*coordinates[m,2])
                end
                
            elseif abs(coordinates[m,2] - coordinates[n,2])==Ny-1 #Magneto Periodic BC
                if coordinates[m,2] > coordinates[n,2]
                    H[m,n] = t*exp(1im*2*pi*α*coordinates[m,1]*Ny)
                elseif coordinates[m,2] < coordinates[n,2]
                    H[m,n] = t*exp(-1im*2*pi*α*coordinates[m,1]*Ny)
                end
                
            # hopping in the bulk
            else 
                if gauge=="Landau"
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = t*exp(1im*2*pi*α*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = t*exp(-1im*2*pi*α*coordinates[m,2])
                    else
                        H[m,n] = t*exp(0)
                    end
                elseif gauge=="Symmetric"
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = t*exp(1im*pi*α*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = t*exp(-1im*pi*α*coordinates[m,2])
                    elseif coordinates[m,2] > coordinates[n,2]
                        H[m,n] = t*exp(-1im*pi*α*coordinates[m,1])
                    elseif coordinates[m,2] < coordinates[n,2]
                        H[m,n] = t*exp(1im*pi*α*coordinates[m,1])
                    else
                        H[m,n] = t*exp(0)
                    end
                end
            end
        end
    end    
    return H
end

function HubbardHofstadter(Lattice, Params, mb_basis)
    
    if HardCore==true
        H_MB = HoppingOp(Lattice, mb_basis, Params)
        H_Total_full = H_MB  
        H_Total_full = (H_Total_full'+H_Total_full)/2
    elseif HardCore==false 
        H_MB = HoppingOp(Lattice, mb_basis, Params)
        H_Int = InteractionOp(pn, Nx, Ny, U, HardCore) 
        H_Total_full = H_MB + H_Int 
        H_Total_full = (H_Total_full'+H_Total_full)/2 
    end
 
    return H_Total_full
end

struct ModelParams
    Nx
    Ny
    α
    gauge
    HardCore
    pn
    U
    perturbation
    imp_str
    Nphi
    Nd
    GroundStateDegeneracy
    lb
    function ModelParams(Nx, Ny, α, pn)
        Nphi = abs(Nx * Ny * α)
        Nd = Int(Nphi-2*pn)
        GroundStateDegeneracy = factorial(Nd + pn - 1) / (factorial(pn - 1)*factorial(Nd)) * (Nphi / pn)
        lb = 1/sqrt(2*pi*abs(α))
        new(Nx, Ny, α, gauge, HardCore, pn, U, perturbation, imp_str, Nphi, Nd, GroundStateDegeneracy, lb)
    end
end