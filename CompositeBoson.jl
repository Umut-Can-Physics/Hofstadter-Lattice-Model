include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots

pn = 2
Nx = 6
Ny = 6
alpha = -1/6
U = 1.0
periodicity = true
HardCore = true
Nphi = abs(Nx * Ny * alpha)
Nd = Int(Nphi-2*pn)
GroundStateDegeneracy = factorial(Nd + pn - 1) / (factorial(pn - 1)*factorial(Nd)) * (Nphi / pn)
Nev = Int(GroundStateDegeneracy) + 1
gauge = "Landau" # "Landau", "Symmetric"
imp_str = 0.1
perturbation = false

problem_type = "MB" # 'SP' or 'MB'
method = "Lapack" # "Lapack", "Arpack", "KrylovKit"
E_mb, ψ_mb, H_mb = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_mb[1:10]))

type = "fermion" # HardCore
UpperLimit = 10
shift_amount = 0
lb = 1/sqrt(2*pi*abs(alpha))
mb = MBBasis(pn, Nx, Ny, HardCore)
OccBasis = MBBasis(pn, Nx, Ny, HardCore).occupations
#OccupiedCoordinates = [ParCoord(Nx, Ny, OccBasis, i, type) for i in eachindex(OccBasis)]
#SiteCoords = SiteCoordinates(Nx, Ny)
#ψ0, ψ1, ψ2  = CompositeBosonMBPart(OccBasis, Nx, Ny, UpperLimit, type)

problem_type = "SP"
alpha_prime = -1/18
#sp_basis = fermionstates(NLevelBasis(Nx*Ny), 1) # I suppose that this is the single particle basis consist of the single particle matrix
E_sp, ψ_sp = Solve(pn, Nx, Ny, alpha_prime, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_sp))

#- LATEST 8/30/25 -#
lb_prime = 1/sqrt(2*pi*abs(alpha_prime))
@showprogress for x1 in 1:Nx
    for x2 in 1:Nx
        for x3 in 1:Nx
            for y1 in 1:Ny
                for y3 in 1:Ny
                    for y2 in 1:Ny
                        z_nu0 = x1 + im*y1
                        z_nu1 = x2 + im*y2
                        z_nu2 = x3 + im*y3
                        Ψ_CBB = CompositeBoson(OccBasis, Nx, Ny, lb, lb_prime, UpperLimit, type, ψ_sp, z_nu0, z_nu1, z_nu2)
                        # Normalize rows
                        Ψ_CBB = normalize.(eachcol(Ψ_CBB))
                        Ψ_ED = ψ_mb[:,1:9]
                        Overlapss = []
                        for k in 1:9
                            for l in 1:9
                                push!(Overlapss,abs(Ψ_CBB[k]'*Ψ_ED[:,l])^2)
                            end
                        end
                        overlaps = reshape(Overlapss, (9,9))
                        HilbertSchmidtNorm(overlaps)
                        println(maximum(overlaps))
                    end
                end
            end
        end
    end
end
#- END -#

# Eigenvalues of the Hamiltonian in the CB basis #

function gram_schmidt_qr(vectors::Vector{Vector{T}}) where T<:Number
    # Convert array-of-vectors into a matrix (columns are the vectors)
    M = hcat(vectors...)  
    
    # QR factorization
    F = qr(M)
    
    # Extract Q (orthonormal basis as a matrix)
    Q = Matrix(F.Q)
    
    # Return as array-of-vectors
    return [Q[:,i] for i in 1:size(Q,2)]
end
Ψ_CBB_orthonormal = gram_schmidt_qr(Ψ_CBB)
H = HubbardHofstadter(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, perturbation, imp_str)
H_new = zeros(ComplexF64, length(Ψ_CBB_orthonormal), length(Ψ_CBB_orthonormal))
N = Nx*Ny
Number_new = zeros(ComplexF64, length(Ψ_CBB_orthonormal), length(Ψ_CBB_orthonormal), N)
for i in eachindex(Ψ_CBB_orthonormal)
    println("Processing state ", i, " out of ", length(Ψ_CBB_orthonormal))
    for j in eachindex(Ψ_CBB_orthonormal)
        H_new[i, j] = Ψ_CBB_orthonormal[i]'*Matrix(H.data)*Ψ_CBB_orthonormal[j]
        for m in 1:N
            Number_new[i, j, m] = Ψ_CBB_orthonormal[i]'*Matrix(number(mb,m).data)*Ψ_CBB_orthonormal[j]
        end
    end
end
E_new, ψ_new = eigen(H_new)
E_new, ψ_new = SortStates(E_new, ψ_new)
println("Eigenvalues in the Composite Boson Basis: ", E_new)
scatter(real(E_new), label="Eigenvalues in the Composite Boson Basis")
scatter!(real(E_mb[1:9]), label="MB Eigenvalues from ED")
Density_new = zeros(ComplexF64,Nx, Ny)
k = 1
for m in 1:N
    Density_new[m] = ψ_new[:,k]'*Number_new[:,:,m]*ψ_new[:,k]
end
heatmap(real(Density_new))

N = Nx*Ny
Density = zeros(N)
ψ'*Number_new.data*ψ

Overlaps_ortonormal = []
for k in 1:9
    for l in 1:9
        push!(Overlaps_ortonormal,abs(Ψ_CBB_orthonormal[k]'*Ψ_ED[:,l])^2)
    end
end
overlaps = reshape(Overlapss, (9,9))


#= ϕ_1 = ψ_sp[:,1]
ϕ_2 = ψ_sp[:,2]

ψ_CB_11 = zeros(ComplexF64, length(OccBasis))
ψ_CB_22 = zeros(ComplexF64, length(OccBasis))
ψ_CB_12 = zeros(ComplexF64, length(OccBasis))
TEST = []
for i in eachindex(OccBasis)
    # coordinates of the occupied sites in the MB basis i
    r = OccupiedCoordinates[i]
    # find the index of these coordinates in the site coordinates
    rr = indexin(r, SiteCoords)

    # corresponding sp basis amplitude
    A1 = ϕ_1[rr]
    A2 = ϕ_2[rr]
    # Composite boson guess
    ψ_CB_11[i] = A1[1] * A1[2] # \phi_1(r_1) \phi_1(r_2) 
    ψ_CB_22[i] = A2[1] * A2[2] # \phi_2(r_1) \phi_2(r_2) 
    ψ_CB_12[i] = A1[1] * A2[2] # \phi_1(r_1) \phi_2(r_2) 

    push!(TEST,[ψ_CB_11[i]*ψ_CB_22[i]*ψ_CB_12[i]])
end
Vecs = [ψ0, ψ1, ψ2]
Ψ_CB = []
for j in 1:3
    push!(Ψ_CB, normalize(ψ_CB_11.*Vecs[j]))
    push!(Ψ_CB, normalize(ψ_CB_22.*Vecs[j]))
    push!(Ψ_CB, normalize(ψ_CB_12.*Vecs[j]))
end

scatter(abs.(Ψ_CB[9]))
scatter(abs.(ψ_mb[:,9]))

heatmap(RealSpaceDensity(Nx, Ny, ψ_mb[:,1], mb))
densityY = zeros(Nx, Ny)
for i in 1:9
    densityY += RealSpaceDensity(Nx, Ny, Ψ_CB[i], mb)
end
heatmap(densityY./9)
sum(densityY)

Ψ_ED = ψ_mb[:,1:9]
Overlapss = []
for k in 1:9
    for l in 1:9
        push!(Overlapss,abs(Ψ_CB[k]'*Ψ_ED[:,l])^2)
    end
end
overlaps = reshape(Overlapss, (9,9))
maximum(overlaps)
heatmap(overlaps, xlabel="Composite Boson States", ylabel="Exact Diagonalization States", title="Overlaps between Composite Boson and Exact Diagonalization States", c=:blues)

HilbertSchmidtNorm(overlaps)


scatter(abs.(Ψ_CB[3]))
scatter!(abs.(Ψ_CBB[3])) =#