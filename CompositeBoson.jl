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
function random_complex_triples(Nx::Int, Ny::Int, num_samples::Int)
    triples = []

    while length(triples) < num_samples
        # generate 3 random complex numbers
        z = Complex{Int}[]
        while length(z) < 3
            x, y = rand(0:Nx-1), rand(0:Ny-1)
            c = complex(x, y)
            if !(c in z) # ensure uniqueness
                push!(z, c)
            end
        end
        push!(triples, tuple(z...))
    end

    return triples
end
samples = random_complex_triples(Nx::Int, Ny::Int, 50)
samples_dict = Dict{Tuple, Float64}()
@showprogress for (id,sample) in enumerate(samples)
    z_nu0, z_nu1, z_nu2 = sample
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
    #HilbertSchmidtNorm(overlaps)
    # Some of elements in the samples dict becomes disappear after initialization dictionary, as I used an array as a value.
    # That is why I changed the value type from Array, where maximum(overlaps) and HS norm are stored, to Float64 and only store the Hilbert Schmidt norm.
    samples_dict[(z_nu0, z_nu1, z_nu2)] = maximum(overlaps)
end
sorted_pairs = sort(collect(samples_dict), by = x -> x[2], rev = true)
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
z_nu0, z_nu1, z_nu2 = 3+4*im, 3+1*im, 0+0*im
Ψ_CBB = CompositeBoson(OccBasis, Nx, Ny, lb, lb_prime, UpperLimit, type, ψ_sp, z_nu0, z_nu1, z_nu2)
# Normalize rows
Ψ_CBB = normalize.(eachcol(Ψ_CBB))
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
scatter(real(E_new), label="Eigenvalues in the Composite Boson Basis")
scatter!(real(E_mb[1:9]), label="MB Eigenvalues from ED")
heatmapList = []
for k in 1:9
    Density_new = zeros(ComplexF64,Nx, Ny)
    for m in 1:N
        Density_new[m] = ψ_new[:,k]'*Number_new[:,:,m]*ψ_new[:,k]
        #Density_new[m] = Ψ_CBB_orthonormal[k]'*number(mb, m).data*Ψ_CBB_orthonormal[k]
    end
    push!(heatmapList, heatmap(real(Density_new)))
end
p = plot(heatmapList..., layout=(3,3), size=(500,500))
savefig(p, "Density_in_Composite_Boson_Basis2.png")

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

# Visualisation of the triplets z_nu0, z_nu1, z_nu2 on the lattice 
z1 = 3+4im
z2 = 3+1im
z3 = 0+0im
a = b = 1/2
τ = im*Ny/Nx
θ1 = []
θ2 = []
θ3 = []
for x in 0:Nx-1
    for y in 0:Ny-1
    z_coord = x + y*im
    push!(θ1, CM_New_3(Nx, Ny, UpperLimit, z_coord, z1))
    push!(θ2, CM_New_3(Nx, Ny, UpperLimit, z_coord, z2))
    push!(θ3, CM_New_3(Nx, Ny, UpperLimit, z_coord, z3))
    end
end
heatmap(reshape(abs.(θ2), (Nx, Ny)))

x_nu0, y_nu0 = (0,0)
x_nu1, y_nu1 = (2,0)
x_nu2, y_nu2 = (5,0)
# Actual lattice points is rotated by 90 degree in the clockwise direction
scatter([x_nu0, x_nu1, x_nu2], [y_nu0, y_nu1, y_nu2], xlims=[0,Nx-1], ylims=[0,Ny-1], xlabel="x", ylabel="y", legend=false)
z_nu_row0 = [0+0im, 2+2im, 4+4im] 

x_nu0, y_nu0 = (3,4)
x_nu1, y_nu1 = (3,1)
x_nu2, y_nu2 = (0,0)
scatter!([x_nu0, x_nu1, x_nu2], [y_nu0, y_nu1, y_nu2], xlims=[0,Nx-1], ylims=[0,Ny-1], xlabel="x", ylabel="y", legend=false)
z_nu_row1 = [0+2im, 2+4im, 4+0im]

x_nu0, y_nu0 = (0,4)
x_nu1, y_nu1 = (2,0)
x_nu2, y_nu2 = (4,2)
scatter!([x_nu0, x_nu1, x_nu2], [y_nu0, y_nu1, y_nu2], xlims=[0,Nx-1], ylims=[0,Ny-1], xlabel="x", ylabel="y", legend=false)
z_nu_row2 = [0+4im, 2+0im, 4+2im]

Ψ_CBB = CompositeBoson(OccBasis, Nx, Ny, lb, lb_prime, UpperLimit, type, ψ_sp, z_nu_row0, z_nu_row1, z_nu_row2)
#Ψ_CBB = CompositeBoson(OccBasis, Nx, Ny, lb, lb_prime, UpperLimit, type, ψ_sp, z_nu_row0[1], z_nu_row0[2], z_nu_row0[3])
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