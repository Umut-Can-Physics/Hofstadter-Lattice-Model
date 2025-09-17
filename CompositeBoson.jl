include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots
include("Run.jl")

Nev = Int(param.GroundStateDegeneracy) + 6
method = "Lapack" # "Lapack", "Arpack", "KrylovKit"

# MB

ϵ, ψ = SolveMatrix(HH, Nev, method)
scatter(real(ϵ[1:Nev]))

# SP

H_sp = HoppingOp(lat, lat.sp_basis, param)
ϵ_sp, ψ_sp = SolveMatrix(H_sp, Nev, method)
scatter(real(ϵ_sp))

# CB

type = "fermion" # HardCore Boson
UpperLimit = 10
alpha_prime = -1/18
lb_prime = 1/sqrt(2*pi*abs(alpha_prime))
WF = "CB" # "Laughlin" or "CB"
Ψ_CB = Ansatz(lat, mb_basis, UpperLimit, param.lb, lb_prime, WF, ψ_sp)

# Normalize column
Ψ_CB = normalize.(eachcol(Ψ_CB))
Ψ_ED = ψ[:,1:9]
OverlapsMatrix = zeros(ComplexF64, 9, 9)
for k in 1:9
    for l in 1:9
        OverlapsMatrix[k,l] = abs(Ψ_CB[k]'*Ψ_ED[:,l])^2
    end
end
HilbertSchmidtNorm(OverlapsMatrix)
maximum(real(OverlapsMatrix))

sort(real(vec(OverlapsMatrix)), rev=true)
idx_sorted = sortperm(real(vec(OverlapsMatrix)), rev=true)
inds = CartesianIndices(OverlapsMatrix)
inds[idx_sorted]

# Eigenvalues of the Hamiltonian in the CB basis #
Ψ_CB_orthonormal = gram_schmidt_qr(Ψ_CB)
H = HubbardHofstadter(pn, Nx, Ny, α, periodicity, gauge, HardCore, U, perturbation, imp_str)
H_new = zeros(ComplexF64, length(Ψ_CB_orthonormal), length(Ψ_CB_orthonormal))
N = Nx*Ny
Number_new = zeros(ComplexF64, length(Ψ_CB_orthonormal), length(Ψ_CB_orthonormal), N)
for i in eachindex(Ψ_CB_orthonormal)
    println("Processing state ", i, " out of ", length(Ψ_CB_orthonormal))
    for j in eachindex(Ψ_CB_orthonormal)
        H_new[i, j] = Ψ_CB_orthonormal[i]'*Matrix(H.data)*Ψ_CB_orthonormal[j]
#=         for m in 1:N
            Number_new[i, j, m] = Ψ_CB_orthonormal[i]'*Matrix(number(mb,m).data)*Ψ_CB_orthonormal[j]
        end  =#
    end
end
E_new, ψ_new = eigen(H_new)
E_new, ψ_new = SortStates(E_new, ψ_new)
scatter(real(E_new), label="Eigenvalues in the Composite Boson Basis", legend=false)
scatter!(real(E_mb), label="MB Eigenvalues from ED")

for k in 1:9
    Density_new = zeros(ComplexF64,Nx, Ny)
    for m in 1:N
        Density_new[m] = ψ_new[:,k]'*Number_new[:,:,m]*ψ_new[:,k]
        #Density_new[m] = Ψ_CBB_orthonormal[k]'*number(mb, m).data*Ψ_CBB_orthonormal[k]
    end
    push!(heatmapList, heatmap(real(Density_new)))
end