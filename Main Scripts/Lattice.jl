using OffsetArrays

function square_lattice(Nx::Int, Ny::Int)    
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

function flipped_rowmajor_reshape(A::AbstractVector, dims::Tuple{Vararg{Int}})
    N = prod(dims)
    if length(A) != N
        throw(ArgumentError("Length of array ($(length(A))) must match product of dimensions ($N)."))
    end

    # Step 1: row-major reshape
    rowmajor = reshape(collect(A), reverse(dims)...) |> permutedims

    # Step 2: flip vertically
    return rowmajor[end:-1:1, :]
end

# Convenience method
flipped_rowmajor_reshape(A::AbstractVector, dims::Int...) = flipped_rowmajor_reshape(A, dims)

function square_lattice_2(Nx::Int, Ny::Int)    
    site_idx = range(1,Nx*Ny) 
    lattice = OffsetArray(flipped_rowmajor_reshape(site_idx, Nx, Ny), 0:Nx-1, 0:Ny-1)
    #lattice = OffsetArray(reshape(site_idx, (Nx,Ny)), 0:Nx-1, 0:Ny-1) |> transpose
    coordinates = []
    for x in 0:Nx-1
        for y in 0:Ny-1
            coordinates = [coordinates; x; y]
        end
    end
    coordinates = reshape(coordinates, (2, Nx*Ny)) |> transpose
    
    return lattice, coordinates
end

function neighbors(Nx::Int, Ny::Int, periodicity::Bool)
    
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

function neighbors_2(Nx::Int, Ny::Int, periodicity::Bool)

    warning("This function is deprecated. Use `neighbors` instead.")
    
    lattice = square_lattice_2(Nx,Ny)[1]
    Neighbors = []

    # Periodicity On
    if periodicity == true
           
        for i in 0:Nx
            for j in 0:Ny
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