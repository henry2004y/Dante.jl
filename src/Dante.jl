module Dante
# Reimplementation of Dante in Julia.

export solve, set_init_Riemann, EulerExact

using Printf, TimerOutputs

# Stimer()re main timer for global timing of functions
const main_timer = TimerOutput()

# Always call timer() timer() hide implementation details
timer() = main_timer

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
function solve(paramFile="run/PARAM.toml")

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
