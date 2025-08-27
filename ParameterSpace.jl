using StatsBase # for gcd

"""
    generate_parameter_space(N, Nphi; maxNxy=20)

Generate all possible parameter tuples (Nx, Ny, p, q) 
satisfying Nphi = Nx * Ny * (p/q), with p and q coprime.

Arguments:
- N      : Total number of particles (currently not used in constraints, but stored)
- Nphi   : Number of flux quanta
- maxNxy : Maximum search bound for Nx and Ny (default: 20)

Returns:
- Vector of tuples (Nx, Ny, p, q), sorted by phi = p/q
"""
function generate_parameter_space(N, Nphi; maxNxy=20)
    results = []

    for Nx in 1:maxNxy
        for Ny in 1:maxNxy
            denom = Nx * Ny
            num   = Nphi
            g     = gcd(num, denom)
            p     = num ÷ g
            q     = denom ÷ g
            if p > 0 && q > 0
                push!(results, (Nx, Ny, p, q))
            end
        end
    end

    # remove duplicates and sort by phi = p/q
    results = unique(results)
    return sort(results, by = x -> x[3] / x[4])
end

N = 2
Nphi = 2
params = generate_parameter_space(N, Nphi; maxNxy=10)

for (Nx, Ny, p, q) in params
    println("Nx=$Nx, Ny=$Ny, Aspect Ratio=$(Nx/Ny) ,p=$p, q=$q, phi=$(p//q): Hilbert Space Dimension=$(binomial(Nx*Ny, N))")
end
