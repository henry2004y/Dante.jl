module Dante
# Reimplementation of finite volume MHD in Julia.

export solve, setInitRiemann, EulerExact

using Printf, TimerOutputs

# Set main timer for global timing of functions
const main_timer = TimerOutput()

# Always call timer() to hide implementation details
timer() = main_timer

include("utility.jl")
include("input.jl")
include("advance.jl")
include("divergence.jl")
include("facevalue.jl")
include("flux.jl")
include("state.jl")
include("boundary.jl")
include("source.jl")
include("time.jl")
include("output.jl")


"""
	solve(paramFile)

Run the model given input parameter file `paramFile`.
"""
function solve(paramFile="examples/PARAM_sods.toml")

	reset_timer!(timer())

	@timeit timer() "input" param = setParameters(paramFile)

	@timeit timer() "init" state_GV = set_init(param)

	@timeit timer() "advance" advance!(param, state_GV)

	# Compare with analytical solution, if possible
	if param.DoPlot && param.IC == "Riemann"
		@timeit timer() "analytic sol" plot_Riemann_exact(param)
	end

	show(timer())
	println()

	# Check if running test
	if occursin("test", paramFile)
		return param, state_GV
	end
end

end
