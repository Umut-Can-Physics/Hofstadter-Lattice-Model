using QuantumOptics

function MBBasis(pn, Nx, Ny, HardCore)
    N = Nx*Ny
    sp_basis = NLevelBasis(N)
    if HardCore==false
        N_States = bosonstates(sp_basis, pn)
        N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    elseif HardCore==true
        N_States = fermionstates(sp_basis, pn)
        N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    end
    return N_Basis_MB
end 