# I don't like the current implementation for timestepping!
"Explicit time advance."
function advance!(param, state_GV)
	(;verbose) = param

	faceState, faceGradient = init_face_value(param)

	faceFluxLR, faceFlux, speedFlux = init_flux(param)

	source_GV, U, div_G = init_source(param)

	time_G = init_timestep(param)

	t, it = 0.0, 0

	verbose && @info "Finished initialization..."

	if param.TimeAccurate # Advance with time
		if param.Order == 1
			while t < param.tEnd
				# Set boundary conditions
				@timeit timer() "BC" set_cell_boundary!(
					param, state_GV)

				# Calculate face value
				@timeit timer() "face value" calc_face_value!(
					param, state_GV, faceState, faceGradient)

				# Calculate face flux
				@timeit timer() "face flux" calc_face_flux!(
					param, faceState, faceFlux, speedFlux, faceFluxLR)

				# Calculate source
				@timeit timer() "source" calc_source!(
					param, state_GV, source_GV, U, div_G)

				# Calculate time step
				@timeit timer() "update dt" dt = calc_timestep!(
					param, speedFlux, time_G)

				if t + dt > param.tEnd dt = param.tEnd - t end

				# Update state
				@timeit timer() "update state" update_state!(
					param, state_GV, dt, faceFlux, source_GV)

				t  += dt
				it += 1

				if param.DoPlot && mod(it, param.PlotInterval) == 0
					@timeit timer() "plotvar" plotvar(param, it, state_GV)
				end

				verbose && @info "it, t = $it, $t"
			end
		elseif param.Order == 2 # 2nd order method
			state1_GV = similar(state_GV)

			while t < param.tEnd
				# 1st stage of modified timestepping

				# Set boundary conditions
				@timeit timer() "BC" set_cell_boundary!(
					param, state_GV)

				# Calculate face value
				@timeit timer() "face value" calc_face_value!(
					param, state_GV, faceState, faceGradient)

				# Calculate face flux
				@timeit timer() "face flux" calc_face_flux!(
					param, faceState, faceFlux, speedFlux, faceFluxLR)

				# Calculate source
				@timeit timer() "source" calc_source!(
					param, state_GV, source_GV, U, div_G)

				# Calculate time step
				@timeit timer() "update dt" dt = calc_timestep!(
					param, speedFlux, time_G)

				# Update state in the 1st stage
				dt *= 0.5
				state1_GV .= state_GV
				@timeit timer() "update state" update_state!(
					param, state1_GV, dt, faceFlux, source_GV)

				# 2nd stage of modified timestepping

				# Set boundary conditions
				@timeit timer() "BC" set_cell_boundary!(
					param, state1_GV)

				# Calculate face value
				@timeit timer() "face value" calc_face_value!(
					param, state1_GV, faceState, faceGradient)

				# Calculate face flux
				@timeit timer() "face flux" calc_face_flux!(
					param, faceState, faceFlux, speedFlux, faceFluxLR)

				# Calculate source
				@timeit timer() "source" calc_source!(
					param, state1_GV, source_GV, U, div_G)

				# Update state
				dt *= 2.0
				@timeit timer() "update state" update_state!(
					param, state_GV, dt, faceFlux, source_GV)

				if t + dt > param.tEnd dt = param.tEnd - t end

				t  += dt
				it += 1

				if param.DoPlot && mod(it, param.PlotInterval) == 0
					@timeit timer() "plotvar" plotvar(param, it, state_GV)
				end

				verbose && @info "it, t = $it, $t"
			end
		else
			@error "Higher Order schemes not yet implemented!"
		end
	else # Advance with steps
		state1_GV = similar(state_GV)

		for iStep = 1:param.nStep
			if param.Order == 1
				# Set boundary conditions
				@timeit timer() "BC" set_cell_boundary!(
					param, state_GV)

				# Calculate face value
				@timeit timer() "face value" calc_face_value!(
					param, state_GV, faceState, faceGradient)

				# Calculate face flux
				@timeit timer() "face flux" calc_face_flux!(
					param, faceState, faceFlux, speedFlux, faceFluxLR)

				# Calculate source
				@timeit timer() "source" calc_source!(
					param, state_GV, source_GV, U, div_G)

				# Calculate time step
				@timeit timer() "update dt" calc_timestep!(
					param, speedFlux, time_G)

				# Update state
				@timeit timer() "update state" update_state!(
					param,state_GV,time_G,faceFlux,source_GV)

				it += 1

				if param.DoPlot && mod(it, param.PlotInterval) == 0
					@timeit timer() "plotvar" plotvar(param, it, state_GV)
				end

				verbose && @info "it = $it"
			elseif param.Order == 2
				# 1st stage of modified timestepping

				# Set boundary conditions
				set_cell_boundary!(param, state_GV)

				# Calculate face value
				calc_face_value!(param, state_GV, faceState, faceGradient)

				# Calculate face flux
				calc_face_flux!(param, faceState, faceFlux, speedFlux, faceFluxLR)

				# Calculate source
				calc_source!(param, state_GV, source_GV, U, div_G)

				# Calculate time step
				calc_timestep!(param, speedFlux, time_G)

				# Update state in the 1st stage
				time_G .*= 0.5
				state1_GV .= state_GV
				update_state!(param, state1_GV, time_G, faceFlux, source_GV)

				# 2nd stage of modified timestepping

				# Set boundary conditions
				set_cell_boundary!(param, state1_GV)

				# Calculate face value
				calc_face_value!(param, state1_GV, faceState, faceGradient)

				# Calculate face flux
				calc_face_flux!(param, faceState, faceFlux, speedFlux, faceFluxLR)

				# Calculate source
				calc_source!(param, state1_GV, source_GV, U, div_G)

				# Update state
				time_G *= 2.0
				update_state!(param, state_GV, time_G, faceFlux, source_GV)

				it += 1

				if param.DoPlot && mod(it, param.PlotInterval) == 0
					plotvar(param, it, state_GV)
				end

				verbose && @info "it = $it"
			end
		end
	end
	return
end
