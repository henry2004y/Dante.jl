module IO

include("EulerExact.jl")

export plotvar, plot_Riemann_exact

using ..Parameters: Param, Param3D, Rho_, Ux_, P_
using ..State: set_init_Riemann
using PyPlot

"""Plot 1D variables along a line."""
function plotvar(param::Param, it::Int, state_GV::Array{Float64,4})

   # Now this only works for 1D x!
   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   plotvar = param.PlotVar
   nG = param.nG

   if plotvar == "rho"
      var = @view state_GV[:,:,:,Rho_]
   elseif plotvar == "ux"
      @views var = state_GV[:,:,:,Ux_]./state_GV[:,:,:,Rho_]
   elseif plotvar == "p"
      var = @view state_GV[:,:,:,P_]
   else
      error("unknown plotting varname!")
   end

   x = param.x[1+nG:end-nG]

   var = @view var[iMin:iMax,jMin:jMax,kMin:kMax] # Remove ghost cells
   var = dropdims(var; dims=(2,3))

   plot(x, var, marker=".")

   title("iStep=$(it)")
   legend(labels=[plotvar])
end

"""Plot the analytical solution of shock tube problem."""
function plot_Riemann_exact(param::Param)
   # Obtain the initial states
   Rho, U, P, tEnd = set_init_Riemann(param.RiemannProblemType)
   # Exact solution
   xe, re, ue, pe, ee, te, Me, se =
   EulerExact(Rho[1], U[1], P[1], Rho[end], U[end], P[end], tEnd, 3)

   # Ee = @. pe/((Parameters.γ-1)*re) + 0.5*ue^2

   if param.PlotVar == "rho"
      plot(xe, re, "--")
   elseif param.PlotVar == "p"
      plot(xe, pe, "--")
   elseif param.PlotVar == "ux"
      plot(xe, ue, "--")
   end

   return
end

end
