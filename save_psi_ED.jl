using MAT

# Save the first two columns of ψ as a MATLAB .mat file
matwrite("psi_ED.mat", Dict("psi_ED" => ψ[:,1:2]))
