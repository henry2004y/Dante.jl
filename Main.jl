# Reimplementation of Dante in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu 08/08/2019

include("Parameters.jl")
include("Divergence.jl")
include("FaceValue.jl")
include("Flux.jl")
include("State.jl")
include("Boundary.jl")
include("Source.jl")
include("Time.jl")

using Printf
using .Parameters, .Divergence, .State, .Boundary, .FaceValue, .Flux, .Source
using .Time

function main()

   param = setParameters()

   state_GV = set_init(param)

   faceState, faceGradient = init_face_value(param)

   faceFluxLR, faceFlux, speedFlux = init_flux(param)

   source_GV = init_source(param)

   time_G = init_timestep(param)

   t, it = 0.0, 0

   println("Finished initialization...")

   if param.TimeAccurate # Advance with time
      if param.Order == 1
         while t < param.tEnd
            # Set boundary conditions
            set_cell_boundary!(param, state_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state_GV,faceState,faceGradient)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            calc_source!(param, state_GV, source_GV)

            # Calculate time step
            dt = calc_timestep!(param, speedFlux, time_G)

            if t + dt > param.tEnd dt = param.tEnd - t end

            # Update state
            update_state!(param,state_GV,dt,faceFlux,source_GV)

            t  += dt
            it += 1

            if mod(it,param.PlotInterval) == 0
               plotvar(param,it,state_GV)
            end

            @printf("it,t=%d, %7.4f\n", it, t)
         end
      elseif param.Order == 2 # 2nd order method
         while t < param.tEnd
            # 1st stage of modified timestepping

            # Set boundary conditions
            set_cell_boundary!(param, state_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state_GV,faceState,faceGradient)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            calc_source!(param, state_GV, source_GV)

            # Calculate time step
            dt = calc_timestep!(param, speedFlux, time_G)

            # Update state in the 1st stage
            dt *= 0.5
            state1_GV = copy(state_GV)
            update_state!(param,state1_GV,dt,faceFlux,source_GV)

            # 2nd stage of modified timestepping

            # Set boundary conditions
            set_cell_boundary!(param, state1_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state1_GV,faceState,faceGradient)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            calc_source!(param, state1_GV, source_GV)

            # Update state
            dt *= 2.0
            update_state!(param,state_GV,dt,faceFlux,source_GV)

            if t + dt > param.tEnd dt = param.tEnd - t end

            t  += dt
            it += 1

            if mod(it,param.PlotInterval) == 0
               plotvar(param,it,state_GV)
            end

            @printf("it,t=%d, %7.4f\n", it, t)
         end
      else
         error("Higher Order schemes not yet implemented!")
      end
   else # Advance with steps
      for iStep = 1:param.nStep
         if param.Order == 1
            # Set boundary conditions
            set_cell_boundary!(param, state_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state_GV,faceState,faceGradient)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            calc_source!(param, state_GV, source_GV)

            # Calculate time step
            calc_timestep!(param, speedFlux, time_G)

            # Update state
            update_state!(param,state_GV,time_G,faceFlux,source_GV)

            it += 1

            if mod(it,param.PlotInterval) == 0
               plotvar(param,it,state_GV)
            end

            @printf("it=%d\n", it)
         elseif param.Order == 2
            # 1st stage of modified timestepping

            # Set boundary conditions
            set_cell_boundary!(param, state_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state_GV,faceState,faceGradient)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            calc_source!(param, state_GV, source_GV)

            # Calculate time step
            calc_timestep!(param, speedFlux, time_G)

            # Update state in the 1st stage
            time_G .*= 0.5
            state1_GV = copy(state_GV)
            update_state!(param,state1_GV,time_G,faceFlux,source_GV)

            # 2nd stage of modified timestepping

            # Set boundary conditions
            set_cell_boundary!(param, state1_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state1_GV,faceState,faceGradient)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            calc_source!(param, state1_GV, source_GV)

            # Update state
            time_G *= 2.0
            update_state!(param,state_GV,time_G,faceFlux,source_GV)

            it += 1

            if mod(it,param.PlotInterval) == 0
               plotvar(param,it,state_GV)
            end

            @printf("it=%d\n", it)
         end
      end
   end

end
