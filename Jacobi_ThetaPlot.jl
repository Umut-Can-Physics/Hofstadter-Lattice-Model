# In Progress
# To check behaviour of the Jacobi Theta Function with arbitrary characteristics

pn = 2
Nx = Ny = 10
HardCore = true
OccBasis = MBBasis(pn, Nx, Ny, HardCore).occupations

a = 1/2
b = -1/2
z = Z(OccBasis, 1, type, shift_amount) / Nx
τ = im*Ny/Nx
    UpperLimit = 1
    ThetaFun = [exp(im*pi*τ*(n+a).^2 + 2*pi*im*(n+a)*(z+b)) for n in -UpperLimit:UpperLimit] 
    scatter(real(ThetaFun))