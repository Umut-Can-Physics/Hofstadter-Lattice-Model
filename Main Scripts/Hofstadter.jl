module Hofstadter

using Revise

# MB Functions
includet("ED.jl")
includet("Lattice.jl")
includet("Model.jl")
includet("Operators.jl")
includet("Utilities.jl")

# Laughlin and CB Functions
includet("../Laughlin Scripts/GeneralizedLaughlin.jl")
includet("../Laughlin Scripts/JacobiThetaFunction.jl")

end