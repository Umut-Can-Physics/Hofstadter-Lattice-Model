function SingleParticleModel(Nx::Int, Ny::Int, α::Float64 ,periodicity::Bool, gauge::String)
    
    if gauge != "Landau" && gauge != "Symmetric"
        error("Invalid gauge choice. Use 'Landau' or 'Symmetric'.")
    end

    @warn("Symmetric gauge is not implemented yet, use Landau gauge instead.")

    neig = neighbors(Nx, Ny, periodicity)
    coordinates = square_lattice(Nx, Ny)[2]
     
    N = Nx*Ny
    t = -1
    H = zeros(Complex{Float64},N,N)
    
    for m in 1:N
        for n in 1:N
            if m in neig[n] 

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
    end
    
    return H
end

function HubbardHofstadter(pn::Int, Nx::Int, Ny::Int, α::Float64, periodicity::Bool, gauge::String, HardCore::Bool, U::Float64, perturbation::Bool, imp_str::Float64)

    if HardCore==true
        H_MB = MBOp(pn, Nx, Ny, α, periodicity, gauge, HardCore, perturbation, imp_str)
        H_Total_full = H_MB  
        H_Total_full = (H_Total_full'+H_Total_full)/2
    elseif HardCore==false
        H_MB = MBOp(pn, Nx, Ny, α, periodicity, gauge, HardCore, perturbation, imp_str)
        H_Int = InteractionOp(pn, Nx, Ny, U, HardCore)
        H_Total_full = H_MB + H_Int
        H_Total_full = (H_Total_full'+H_Total_full)/2
    end

    return H_Total_full
end