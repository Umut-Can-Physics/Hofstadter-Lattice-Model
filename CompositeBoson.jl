include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots

pn = 2
Nx = 10
Ny = 10
alpha = -3/50
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
method = "Arpack" # "Lapack", "Arpack", "KrylovKit"
E_mb, ψ_mb = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_mb))

type = "fermion" # HardCore
UpperLimit = 10
shift_amount = 0
lb = 1/sqrt(2*pi*abs(alpha))
mb = MBBasis(pn, Nx, Ny, HardCore)
OccBasis = MBBasis(pn, Nx, Ny, HardCore).occupations
OccupiedCoordinates = [ParCoord(Nx, Ny, OccBasis, i, type) for i in eachindex(OccBasis)]
SiteCoords = SiteCoordinates(Nx, Ny)
ψ0, ψ1, ψ2  = CompositeBosonMBPart(OccBasis, Nx, Ny, UpperLimit, type)

problem_type = "SP"
alpha = -1/50 
lb_prime = 1/sqrt(2*pi*abs(alpha))
#sp_basis = fermionstates(NLevelBasis(Nx*Ny), 1) # I suppose that this is the single particle basis consist of the single particle matrix
E_sp, ψ_sp = Solve(pn, Nx, Ny, alpha, periodicity, gauge, HardCore, U, Nev, perturbation, imp_str, method, problem_type)
scatter(real(E_sp))

#- LATEST 8/30/25 -#
Ψ_CBB, TESTT = CompositeBoson(OccBasis, Nx, Ny, lb, lb_prime, UpperLimit, type, ψ_sp)
# Normalize rows
Ψ_CBB = normalize.(eachcol(Ψ_CBB))
scatter(abs.(Ψ_CBB[1]))
Ψ_ED = ψ_mb[:,1:9]
Overlapss = []
for k in 1:9
    for l in 1:9
        push!(Overlapss,abs(Ψ_CBB[k]'*Ψ_ED[:,l])^2)
    end
end
overlaps = reshape(Overlapss, (9,9))
maximum(overlaps)
heatmap(overlaps, xlabel="Composite Boson States", ylabel="Exact Diagonalization States", title="Overlaps between Composite Boson and Exact Diagonalization States", c=:blues)
#- END -#


ϕ_1 = ψ_sp[:,1]
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
scatter!(abs.(Ψ_CBB[3]))