using OffsetArrays

function square_lattice(Nx,Ny)    
    site_idx = range(1,Nx*Ny) 
    lattice = OffsetArray(reshape(site_idx, (Nx,Ny)), 0:Nx-1, 0:Ny-1) |> transpose
    coordinates = []
    for y in 0:Ny-1
        for x in 0:Nx-1
            coordinates = [coordinates; x; y]
        end
    end
    coordinates = reshape(coordinates, (2, Nx*Ny)) |> transpose
    
    return lattice, coordinates
end

function neighbors(Nx, Ny, periodicity)
    
    lattice = square_lattice(Nx,Ny)[1]
    Neighbors = []

    # Periodicity On
    if periodicity == true
           
        for j in 0:Ny-1
            for i in 0:Nx-1
                x = [lattice[mod(j,Ny),mod(i-1,Nx)],lattice[mod(j+1,Ny),mod(i,Nx)],lattice[mod(j,Ny),mod(i+1,Nx)],lattice[mod(j-1,Ny),mod(i,Nx)]]
                x = unique(x)
                push!(Neighbors,x)
            end
            
        end
    # Periodicity Off (Hard-Wall)
    elseif periodicity == false
        
        for j in 0:Ny-1
            for i in 0:Nx-1
                if j == 0 || i == 0 || j == Ny-1 || i == Nx-1 
                    new_neighbors = []
                    if j != 0
                        push!(new_neighbors, lattice[j-1,i])  
                    end
                    if i != 0
                        push!(new_neighbors, lattice[j,i-1])  
                    end
                    if j != Ny-1
                        push!(new_neighbors, lattice[j+1,i])  
                    end
                    if i != Nx-1
                        push!(new_neighbors, lattice[j,i+1])  
                    end
                else
                    new_neighbors = [
                        lattice[j,i-1],
                        lattice[j+1,i],
                        lattice[j,i+1],
                        lattice[j-1,i]
                        ]
                    push!(Neighbors,new_neighbors)
                end
            Neighbors = push!(Neighbors,new_neighbors)
            Neighbors = unique(Neighbors)
            end
        end
        
    end
    
    return Neighbors
end