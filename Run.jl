include("Main Scripts/Hofstadter.jl")
using .Hofstadter
using Plots

Nx = 4
Ny = 4
periodicity = true
lat = Lattice(Nx, Ny, periodicity)

pn = 2
α = -1/4
U = 0
HardCore = true
imp_str = 0.1
perturbation = false
gauge = "Landau" # "Landau", "Symmetric"
param = ModelParams(Nx, Ny, α, pn)

mb_basis = MBBasis(lat.sp_basis, param)

HH = HubbardHofstadter(lat, param, mb_basis)