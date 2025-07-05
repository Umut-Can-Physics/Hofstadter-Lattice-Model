function SingleParticleModel(Nx, Ny, alpha ,periodicity)
    
    neig = neighbors(Nx, Ny, periodicity)
    coordinates = square_lattice(Nx, Ny)[2]
    
    N = Nx*Ny
    t = -1
    H = zeros(Complex{Float64},Nx*Ny,Nx*Ny)
    
    for m in 1:N
        for n in 1:N
            if m in neig[n] 

                if abs(coordinates[m,1]-coordinates[n,1])==Nx-1
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,2])
                    end
                    
                elseif abs(coordinates[m,2]-coordinates[n,2])==Ny-1 #Magneto Periodic BC
                    if coordinates[m,2] > coordinates[n,2]
                        H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,1]*Ny)
                    elseif coordinates[m,2] < coordinates[n,2]
                        H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,1]*Ny)
                    end
                    
                else
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,2])
                    else
                        H[m,n] = t*exp(0)
                    end
                    
                end
            else
                
                H[m,n] = 0
            end
        end
    end
    
    return H
end

function HofstadterHubbard(pn, Nx, Ny, alpha, periodicity, HardCore, U, perturbation, imp_str)

    if HardCore==true
        H_MB = MBOp(pn, Nx, Ny, alpha, periodicity, HardCore, perturbation, imp_str)
        H_Total_full = H_MB 
        H_Total_full = (H_Total_full'+H_Total_full)/2
        
    elseif HardCore==false
        H_MB = MBOp(pn, Nx, Ny, alpha, periodicity, HardCore, perturbation, imp_str)
        H_Int = InteractionOp(pn, Nx, Ny, U, HardCore)
        H_Total_full = H_MB + H_Int
        H_Total_full = (H_Total_full'+H_Total_full)/2
    end

    return H_Total_full
end