module Dante
# Reimplementation of Dante in Julia.

export main, set_init_Riemann, EulerExact

using Printf, LinearAlgebra

include("Parameters.jl")
include("Advance.jl")
include("Divergence.jl")
include("FaceValue.jl")
include("Flux.jl")
include("State.jl")
include("Boundary.jl")
include("Source.jl")
include("Time.jl")
include("IO.jl")

using .Parameters, .Divergence, .State, .Boundary, .FaceValue, .Flux, .Source,
	.Time, .IO

"""
	main(paramFile)

Main file for execution.
"""
function main(paramFile="run/PARAM.toml")

	param = setParameters(paramFile)

	state_GV = set_init(param)

	advance!(param, state_GV)

	# Compare with analytical solution, if possible
	if param.DoPlot && param.IC == "Riemann"
		plot_Riemann_exact(param)
	end

	# Check if running test
	if occursin("test", paramFile)
		return param, state_GV
	end

end

end
