include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots

pn = 2
Nx = 6
Ny = 6
α = -1/6
U = 1.0
periodicity = true
HardCore = true
Nphi = abs(Nx * Ny * α)
Nd = Int(Nphi-2*pn)
GroundStateDegeneracy = factorial(Nd + pn - 1) / (factorial(pn - 1)*factorial(Nd)) * (Nphi / pn)
Nev = Int(GroundStateDegeneracy) + 1
gauge = "Landau" # "Landau", "Symmetric"
imp_str = 0.1
perturbation = false

problem_type = "MB" # 'SP' or 'MB'
method = "Arpack" # "Lapack", "Arpack", "KrylovKit"
E_mb, ψ_mb, H_mb = Solve(pn, Nx, Ny, α, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_mb[1:Nev]))

type = "fermion" # HardCore Boson
UpperLimit = 10
shift_amount = 0
lb = 1/sqrt(2*pi*abs(α))
mb = MBBasis(pn, Nx, Ny, HardCore)
OccBasis = MBBasis(pn, Nx, Ny, HardCore).occupations

problem_type = "SP"
alpha_prime = -1/18
E_sp, ψ_sp = Solve(pn, Nx, Ny, alpha_prime, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_sp))

lb_prime = 1/sqrt(2*pi*abs(alpha_prime))
WF = "CB" # "Laughlin" or "CB"
Ψ_CB = Ansatz(OccBasis, Nx, Ny, α, lb, lb_prime, UpperLimit, type, shift_amount, WF, ψ_sp)
# Normalize column
Ψ_CB = normalize.(eachcol(Ψ_CB))
Ψ_ED = ψ_mb[:,1:9]
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