include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots
using Serialization
include("Run.jl")

Nev = Int(param.GroundStateDegeneracy) + 10
method = "Arpack" # "Lapack", "Arpack", "KrylovKit"

# MB

ϵ, ψ = SolveMatrix(HH, Nev, method)
scatter(real(ϵ[1:Nev]))

# SP
alpha_prime = -1/18
param_CB = ModelParams(Nx, Ny, alpha_prime, 1)
H_sp = HoppingOp(lat, lat.sp_basis, param_CB)
ϵ_sp, ψ_sp = SolveMatrix(H_sp, Nev, method)
scatter(real(ϵ_sp))

# Analytic Ansatz for SP
UpperLimit = 10
τ = im*Ny/Nx
Nphi_prime = abs(Int(Nx*Ny*alpha_prime))
ψ_sp_1 = SpAnalyticWaveFunction(1, Nx, Ny, lat, UpperLimit, Nphi_prime, τ)
ψ_sp_2 = SpAnalyticWaveFunction(2, Nx, Ny, lat, UpperLimit, Nphi_prime, τ)
ψ_sp_analytic = hcat(ψ_sp_1, ψ_sp_2)

# CB
type = "fermion" # HardCore Boson
WF = "CB" # "Laughlin" or "CB"
qq = 2
pp = nothing
number_of_ansatz = 3*6

Ψ_CB = Ansatz(lat, mb_basis, UpperLimit, param.lb, param_CB.lb, WF, ψ_sp_analytic, qq, number_of_ansatz, pp)

# Normalize column
Ψ_CB = normalize.(eachcol(Ψ_CB))
Ψ_ED = ψ[:,1:number_of_ansatz]
#= OverlapsMatrix = zeros(ComplexF64, number_of_ansatz, number_of_ansatz)
for k in 1:number_of_ansatz
    for l in 1:number_of_ansatz
        OverlapsMatrix[k,l] = Ψ_CB[k]'*Ψ_ED[:,l]
    end
end =#
#HilbertSchmidtNorm(OverlapsMatrix)
#maximum(real(OverlapsMatrix))

# Eigenvalues of the Hamiltonian in the CB basis #
Ψ_CB_orthonormal = gram_schmidt_qr(Ψ_CB)
#Ψ_CB_orthonormal = Ψ_CB
#rank(hcat(Ψ_CB...))
#H = HubbardHofstadter(pn, Nx, Ny, α, periodicity, gauge, HardCore, U, perturbation, imp_str)
H_new = zeros(ComplexF64, length(Ψ_CB_orthonormal), length(Ψ_CB_orthonormal))
N = Nx*Ny
Number_new = zeros(ComplexF64, length(Ψ_CB_orthonormal), length(Ψ_CB_orthonormal), N)
for i in eachindex(Ψ_CB_orthonormal)
    #println("Processing state ", i, " out of ", length(Ψ_CB_orthonormal))
    for j in eachindex(Ψ_CB_orthonormal)
        H_new[i, j] = Ψ_CB_orthonormal[i]'*Matrix(HH.data)*Ψ_CB_orthonormal[j]
#=         for m in 1:N
            Number_new[i, j, m] = Ψ_CB_orthonormal[i]'*Matrix(number(mb,m).data)*Ψ_CB_orthonormal[j]
        end  =#
    end
end 
E_new, ψ_new = eigen(H_new)
E_new, ψ_new = SortStates(E_new, ψ_new)
scatter(real(E_new[1:number_of_ansatz]), label="Eigenvalues in the Composite Boson Basis", legend=false)
scatter!(real(ϵ)[1:10], label="MB Eigenvalues from ED")

#= energy_list = []
for i in 1:50
    push!(energy_list, energy_and_param_list[i][1])
end
minimum(real(energy_list))

for k in 1:9
    Density_new = zeros(ComplexF64,Nx, Ny)
    for m in 1:N
        Density_new[m] = ψ_new[:,k]'*Number_new[:,:,m]*ψ_new[:,k]
        #Density_new[m] = Ψ_CBB_orthonormal[k]'*number(mb, m).data*Ψ_CBB_orthonormal[k]
    end
    push!(heatmapList, heatmap(real(Density_new)))
end =#