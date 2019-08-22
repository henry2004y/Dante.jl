# Reimplementation of Dante in Julia.
#
# This is, to my surprise, not faster than Matlab in some cases. Several
# possible reasons:
# 1. Matlab is somehow parallelized
# 2. Matlab is hightly optimmized in terms of matrix operations
# 3. My Julia implementation is not optimmized.
# Similar behavior is seen in rewriting Kempo1d.
#
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

   faceFluxLR, faceFlux, speedFlux = init_flux(param)

   t, it = 0.0, 0

   println("Finished initialization...")

   if param.TimeAccurate # Advance with time
      if param.Order == 1
         while t < param.tEnd
            # Set boundary conditions
            set_cell_boundary!(param, state_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state_GV)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            source_GV = calc_source(param, state_GV)

            # Calculate time step
            dt = calc_timestep(param, speedFlux)

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
      elseif param.Order == 2 # % 2nd order method
         while t < param.tEnd
            # 1st stage of modified timestepping

            # Set boundary conditions
            set_cell_boundary!(param, state_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state_GV)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            source_GV = calc_source(param, state_GV)

            # Calculate time step
            dt = calc_timestep(param, speedFlux)

            # Update state in the 1st stage
            dt *= 0.5
            state1_GV = copy(state_GV)
            update_state!(param,state1_GV,dt,faceFlux,source_GV)

            # 2nd stage of modified timestepping

            # Set boundary conditions
            set_cell_boundary!(param, state1_GV)

            # Calculate face value
            faceValue = calc_face_value(param, state1_GV)

            # Calculate face flux
            calc_face_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

            # Calculate source
            source_GV = calc_source(param, state1_GV)

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

   end

end
