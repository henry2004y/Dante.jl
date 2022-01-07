# State variables

function set_init(param::Param)
   (;FullSize, nVar, nI, nJ, nK, nG) = param

   state_GV = zeros(FullSize[1],FullSize[2],FullSize[3], nVar)

   density  = @view state_GV[:,:,:,Rho_]
   velocity = @view state_GV[:,:,:,U_]
   B        = @view state_GV[:,:,:,B_]
   pressure = @view state_GV[:,:,:,P_]

   if param.IC == "contact discontinuity"
      density[1:floor(Int,nI/2),:,:] .= 2.0
      density[floor(Int,nI/2+1):end,:,:] .= 1.0

      velocity .= 0.0
      velocity .*= density
      pressure .= 1.0

   elseif param.IC == "density wave"
      density .= 1.0
      density[floor(Int,nI/2):floor(Int,nI/2),:,:] .= 2.0
      velocity[:,:,:,1] .= 1.0
      velocity .*= density
      pressure .= 0.01
   elseif param.IC == "shocktube"

   elseif param.IC == "square wave"
      density .= 1.0
      density[floor(Int,nI/2)-10:floor(Int,nI/2),:,:] .= 2.0
      velocity[:,:,:,1] .= 1.0
      velocity .*= density
      pressure .= 0.01
   elseif param.IC == "Riemann"
      rho, u, p, tEnd, CFL = setInitRiemann(param.RiemannProblemType)
      # Initial Condition for 1D domain
      # Density
      density[1:floor(Int,nI/2)+nG+1,:,:] .= rho[1] # region 1
      density[floor(Int,nI/2)+nG+2:end,:,:] .= rho[2] # region 2
      # Velocity in x
      velocity[1:floor(Int,nI/2)+nG+1,:,:,1] .= u[1] # region 1
      velocity[floor(Int,nI/2)+nG+2:end,:,:,1] .= u[2] # region 2
      # Pressure
      pressure[1:floor(Int,nI/2)+nG+1,:,:] .= p[1] # region 1
      pressure[floor(Int,nI/2)+nG+2:end,:,:] .= p[2] # region 2

   else
      error("unknown initial condition type!")
   end

   return state_GV
end

"""
	setInitRiemann(RiemannProblemType, Verbose)

Set the initial conditions of Riemann problems. Note that currently tEnd and CFL number can
only be set in PARAM.toml, and cannot be changed afterwards.
"""
function setInitRiemann(RiemannProblemType, Verbose=true)

   if RiemannProblemType == 1
      Verbose && println("Case 1: Sods problem")
      p   = [1.0, 0.1]
      u   = [0.0, 0.0]
      rho = [1.0, 0.125]
      tEnd, CFL = 0.1, 0.9 # I need to think of how to set this!!!
   elseif RiemannProblemType == 2
      Verbose && println("Case 2: Left Expansion and right strong shock")
      p   = [1000.0, 0.1]
      u   = [0.0   , 0.0]
      rho = [3.0   , 2.0]
      tEnd, CFL = 0.02, 0.9
   elseif RiemannProblemType == 3
      Verbose && println("Case 3: Right Expansion and left strong shock")
      p   = [7.0, 10.0]
      u   = [0.0, 0.0 ]
      rho = [1.0, 1.0 ]
      tEnd, CFL = 0.1, 0.9
   elseif RiemannProblemType == 4
      Verbose && println("Case 4: Double Shock")
      p   = [450.0, 45.0]
      u   = [20.0 , -6.0]
      rho = [6.0  ,  6.0]
      tEnd, CFL = 0.01, 0.90
   elseif RiemannProblemType == 5
      Verbose && println("Case 5: Double Expansion")
      p   = [40.0,  40.0]
      u   = [-2.0,  2.0 ]
      rho = [1.0 ,  2.5 ]
      tEnd, CFL = 0.03, 0.90
   elseif RiemannProblemType == 6
      Verbose && println("Case 6: Cavitation")
      p   = [0.4 , 0.4]
      u   = [-2.0, 2.0]
      rho = [ 1.0, 1.0]
      tEnd, CFL = 0.1, 0.90
   elseif RiemannProblemType == 7
      Verbose && println("Shocktube problem of G.A. Sod, JCP 27:1, 1978")
      p   = [1.0 , 0.1  ]
      u   = [0.75, 0.0  ]
      rho = [1.0 , 0.125]
      tEnd, CFL = 0.17, 0.90
   elseif RiemannProblemType == 8
      Verbose &&
         println("Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997")
      p   = [3.528, 0.571]
      u   = [0.698, 0    ]
      rho = [0.445, 0.5  ]
      tEnd, CFL = 0.15, 0.90
   elseif RiemannProblemType == 9
      Verbose &&
         println("Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997")
      p   = [10.333,  1.0 ]
      u   = [ 0.92 ,  3.55]
      rho = [3.857 ,  1.0 ]
      tEnd, CFL = 0.09, 0.90
   elseif RiemannProblemType == 10
      Verbose && println("Shocktube problem with supersonic zone")
      p   = [1.0,  0.02]
      u   = [0.0,  0.00]
      rho = [1.0,  0.02]
      tEnd, CFL = 0.162, 0.90
   elseif RiemannProblemType == 11
      Verbose && println("Contact discontinuity")
      p   = [0.5, 0.5]
      u   = [0.0, 0.0]
      rho = [1.0, 0.6]
      tEnd, CFL = 1.0, 0.90
   elseif RiemannProblemType == 12
      Verbose && println("Stationary shock")
      p   = [ 1.0,  0.1 ]
      u   = [-2.0, -2.0 ]
      rho = [ 1.0, 0.125]
      tEnd, CFL = 0.1, 0.28
   else
      error("RiemannProblemType not known!")
   end
   # Print for Riemann Problems
   if Verbose
      println("")
      println("density (L): $(rho[1])")
      println("velocity(L): $(u[1])")
      println("Pressure(L): $(p[1])")
      println("")
      println("density (R): $(rho[2])")
      println("velocity(R): $(u[2])")
      println("Pressure(R): $(p[2])")
   end

   return rho, u, p, tEnd, CFL
end

"Time-accurate state update."
function update_state!(param::Param, state_GV, dt::AbstractFloat, faceFlux::FaceFlux,
   source_GV)

   (;Flux_XV, Flux_YV, Flux_ZV) = faceFlux
   (;CellSize_D, iMin, iMax, jMin, jMax, kMin, kMax, nI, nJ, nK, nG, nVar) = param

   if param.TypeGrid == "Cartesian"
      # No need for volume and face if the grid is uniform Cartesian
      if !param.UseConservative
         @inbounds for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI
            state_GV[i+nG,j+nG,k+nG,iVar] -= dt*( -source_GV[i,j,k,iVar] +
            (Flux_XV[i+1,j,k,iVar] - Flux_XV[i,j,k,iVar])/CellSize_D[1] +
            (Flux_YV[i,j+1,k,iVar] - Flux_YV[i,j,k,iVar])/CellSize_D[2] +
            (Flux_ZV[i,j,k+1,iVar] - Flux_ZV[i,j,k,iVar])/CellSize_D[3])
         end
      else
         @inbounds for k=kMin:kMax, j=jMin:jMax, i=iMin:iMax
            u = □(state_GV[i,j,k,Ux_], state_GV[i,j,k,Uy_], state_GV[i,j,k,Uz_])
            b = □(state_GV[i,j,k,Bx_], state_GV[i,j,k,By_], state_GV[i,j,k,Bz_])

            state_GV[i,j,k,E_] = state_GV[i,j,k,P_]/(γ-1.0) +
               0.5/state_GV[i,j,k,Rho_]*u + 0.5*b
         end

         @inbounds for iVar=Rho_:Bz_, k=1:nK, j=1:nJ, i=1:nI
            state_GV[i+nG,j+nG,k+nG,iVar] -= dt*(source_GV[i,j,k,iVar] +
            (Flux_XV[i+1,j,k,iVar] - Flux_XV[i,j,k,iVar])/CellSize_D[1] +
            (Flux_YV[i,j+1,k,iVar] - Flux_YV[i,j,k,iVar])/CellSize_D[2] +
            (Flux_ZV[i,j,k+1,iVar] - Flux_ZV[i,j,k+1,iVar])/CellSize_D[3])
         end

         @inbounds for k=1:nK, j=1:nJ, i=1:nI
            u = □(state_GV[i+nG,j+nG,k+nG,Ux_],
               state_GV[i+nG,j+nG,k+nG,Uy_],
               state_GV[i+nG,j+nG,k+nG,Uz_])
            b = □(state_GV[i+nG,j+nG,k+nG,Bx_],
               state_GV[i+nG,j+nG,k+nG,By_],
               state_GV[i+nG,j+nG,k+nG,Bz_])

            state_GV[i+nG,j+nG,k+nG,E_] -= dt*( -source_GV[i,j,k,E_] +
               (Flux_XV[i+1,j,k,E_] - Flux_XV[i,j,k,E_]) / CellSize_D[1] +
               (Flux_YV[i,j+1,k,E_] - Flux_YV[i,j,k,E_]) / CellSize_D[2] +
               (Flux_ZV[i,j,k+1,E_] - Flux_ZV[i,j,k,E_]) / CellSize_D[3] )

            state_GV[i+nG,j+nG,k+nG,P_] = (γ-1.0)*(state_GV[i+nG,j+nG,k+nG,E_] -
               0.5/state_GV[i+nG,j+nG,k+nG,Rho_]*u - 0.5*b)
         end
      end
   else
      # Need volume and face
      state_GV .= 0.0
   end
   return
end

"Local timestepping state update."
function update_state!(param::Param, state_GV, dt, faceFlux::FaceFlux, source_GV)
   (;Flux_XV, Flux_YV, Flux_ZV) = faceFlux
   (;CellSize_D, iMin, iMax, jMin, jMax, kMin, kMax, nI, nJ, nK, nG, nVar) = param

   if param.TypeGrid == "Cartesian"
      # No need for volume and face if the grid is uniform Cartesian
      if !param.UseConservative
         @inbounds for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI
            state_GV[i+nG,j+nG,k+nG,iVar] -= dt[i,j,k]*(-source_GV[i,j,k,iVar] +
            (Flux_XV[i+1,j,k,iVar] - Flux_XV[i,j,k,iVar])/CellSize_D[1] +
            (Flux_YV[i,j+1,k,iVar] - Flux_YV[i,j,k,iVar])/CellSize_D[2] +
            (Flux_ZV[i,j,k+1,iVar] - Flux_ZV[i,j,k,iVar])/CellSize_D[3])
         end
      else
         @inbounds for k=kMin:kMax, j=jMin:jMax, i=iMin:iMax
            u = □(state_GV[i,j,k,Ux_], state_GV[i,j,k,Uy_], state_GV[i,j,k,Uz_])
            b = □(state_GV[i,j,k,Bx_], state_GV[i,j,k,By_], state_GV[i,j,k,Bz_])

            state_GV[i,j,k,E_] = state_GV[i,j,k,P_]/(γ-1.0) +
               0.5/state_GV[i,j,k,Rho_]*u + 0.5*b
         end

         @inbounds for iVar=Rho_:Bz_, k=1:nK, j=1:nJ, i=1:nI
            state_GV[i+nG,j+nG,k+nG,iVar] -= dt[i,j,k]*( source_GV[i,j,k,iVar] +
            (Flux_XV[i+1,j,k,iVar] - Flux_XV[i,j,k,iVar]) / CellSize_D[1] +
            (Flux_YV[i,j+1,k,iVar] - Flux_YV[i,j,k,iVar]) / CellSize_D[2] +
            (Flux_ZV[i,j,k+1,iVar] - Flux_ZV[i,j,k+1,iVar]) / CellSize_D[3] )
         end

         @inbounds for k=1:nK, j=1:nJ, i=1:nI
            u = □(state_GV[i+nG,j+nG,k+nG,Ux_],
               state_GV[i+nG,j+nG,k+nG,Uy_],
               state_GV[i+nG,j+nG,k+nG,Uz_])
            b = □(state_GV[i+nG,j+nG,k+nG,Bx_],
               state_GV[i+nG,j+nG,k+nG,By_],
               state_GV[i+nG,j+nG,k+nG,Bz_])

            state_GV[i+nG,j+nG,k+nG,E_] -= dt[i,j,k]*( -source_GV[i,j,k,E_] +
               (Flux_XV[i+1,j,k,E_] - Flux_XV[i,j,k,E_]) / CellSize_D[1] +
               (Flux_YV[i,j+1,k,E_] - Flux_YV[i,j,k,E_]) / CellSize_D[2] +
               (Flux_ZV[i,j,k+1,E_] - Flux_ZV[i,j,k,E_]) / CellSize_D[3] )

            state_GV[i+nG,j+nG,k+nG,P_] = (γ-1.0)*(state_GV[i+nG,j+nG,k+nG,E_] -
               0.5/state_GV[i+nG,j+nG,k+nG,Rho_]*u - 0.5*b)
         end
      end
   else
      # Need volume and face
      state_GV .= 0.0
   end
   return
end