using ProgressMeter
using SparseArrays

# Convert matrices to operator type 
function SPOp(Nx::Int, Ny::Int, α::Float64, periodicity::Bool, gauge::String)
    
    N = Nx*Ny

    sp_basis = NLevelBasis(N) 

    sp_matrix = SingleParticleModel(Nx, Ny, α ,periodicity, gauge)

    H = SparseOperator(sp_basis)

    N, = size(sp_matrix)
    
    for m in 1:N
        for n in 1:N
            H += sp_matrix[m,n] * transition(sp_basis, m, n)
        end
    end
    
    return dense((H'+H)/2)
end

function MBOp(pn::Int, Nx::INt, Ny::Int, α::Float64, periodicity::Bool, gauge::String, HardCore::Bool, perturbation::Bool, imp_str::Float64)

    mb_basis = MBBasis(pn, Nx, Ny, HardCore)

    sp_op = SPOp(Nx, Ny, α, periodicity, gauge)

    mb_op = SparseOperator(mb_basis)
    N = size(sp_op.data, 1) # More robust way to get dimension

    # Pre-allocate lists to build the sparse matrix components directly
    I = Vector{Int}()
    J = Vector{Int}()
    V = Vector{ComplexF64}()

    println("Building many-body operator...")
    @showprogress for j in 1:N
        for i in 1:N
            if !iszero(sp_op.data[i, j])
                transition_op = transition(mb_basis, i, j)
                rows, cols, vals = findnz(transition_op.data)
                scale_factor = sp_op.data[i, j]
                for k in eachindex(vals)
                    push!(I, rows[k])
                    push!(J, cols[k])
                    push!(V, vals[k] * scale_factor)
                end
            end
        end
    end

    # Construct the sparse matrix directly
    sparse_data = sparse(I, J, V, only(mb_basis.shape), only(mb_basis.shape))
    mb_op.data = sparse_data

    if perturbation == true
        mb_op += imp_str*number(mb_basis, 1)
    end

    return mb_op
end

function InteractionOp(pn::Int, Nx::Int, Ny::Int, U::Float64, HardCore::Bool)

    N = Nx*Ny
    sp_basis = NLevelBasis(N)

    mb_basis = MBBasis(pn, Nx, Ny, HardCore)

    basis2 = sp_basis ⊗ sp_basis
    
    Vint2 = SparseOperator(basis2)

    for n in 1:N
        Vint2 += U/2*transition(sp_basis,n,n) ⊗ transition(sp_basis,n,n)
    end

    Vint_mb = manybodyoperator(mb_basis, Vint2)

    return Vint_mb
end