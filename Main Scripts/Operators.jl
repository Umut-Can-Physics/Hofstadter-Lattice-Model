using ProgressMeter
using SparseArrays
using QuantumOptics

# Convert matrices to operator type 
function HoppingOp(Lattice, basis, Params)

    sp_matrix = SingleParticleModel(Lattice, Params.α, Params.gauge)

    H = SparseOperator(basis)
    
    for m in 1:Lattice.N
        for n in Lattice.neig[m]
            H += sp_matrix[m,n] * transition(basis, m, n)
        end
    end
    
    return dense((H'+H)/2)
end

function MBBasis(sp_basis, Params)
    if HardCore==false 
        N_States = bosonstates(sp_basis, Params.pn)
        N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    elseif HardCore==true
        N_States = fermionstates(sp_basis, Params.pn)
        N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    end
    return N_Basis_MB
end 

function InteractionOp(pn::Int, Nx::Int, Ny::Int, U::Float64, HardCore::Bool)
 
    N = Nx*Ny
    sp_basis = NLevelBasis(N)

    mb_basis = MBBasis(pn, N, HardCore)

    basis2 = sp_basis ⊗ sp_basis
    
    Vint2 = SparseOperator(basis2)

    for n in 1:N
        Vint2 += U/2*transition(sp_basis,n,n) ⊗ transition(sp_basis,n,n)
    end

    Vint_mb = manybodyoperator(mb_basis, Vint2)

    return Vint_mb
end