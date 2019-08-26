module Source

export calc_source!, init_source

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_
using ..Parameters: γ
using ..Divergence: divergence!

function init_source(param::Param)
   GridSize = param.GridSize
   nVar = param.nVar
   source_GV = Array{Float64,4}(undef, GridSize..., nVar)

   return source_GV
end

function calc_source!(param::Param, state_GV::Array{Float64,4},
   source_GV::Array{Float64,4})

   nVar, nG = param.nVar, param.nG
   x, y, z = param.x, param.y, param.z

   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   # I can preallocate div_G and DivU!
   nI,nJ,nK = param.nI, param.nJ, param.nK

   div_G = Array{Float64,3}(undef, nI, nJ, nK)
   # Calculate divergence of B using central difference
   B = @view state_GV[:,:,:,B_]
   divergence!(param, B, div_G)

   source_GV[:,:,:,Rho_] .= 0.0

   @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI
      source_GV[i,j,k,Ux_] = -state_GV[i+nG,j+nG,k+nG,Bx_]*div_G[i,j,k]
      source_GV[i,j,k,Uy_] = -state_GV[i+nG,j+nG,k+nG,By_]*div_G[i,j,k]
      source_GV[i,j,k,Uz_] = -state_GV[i+nG,j+nG,k+nG,Bz_]*div_G[i,j,k]
      source_GV[i,j,k,Bx_] =  state_GV[i+nG,j+nG,k+nG,Ux_]*div_G[i,j,k]
      source_GV[i,j,k,By_] =  state_GV[i+nG,j+nG,k+nG,Uy_]*div_G[i,j,k]
      source_GV[i,j,k,Bz_] =  state_GV[i+nG,j+nG,k+nG,Uz_]*div_G[i,j,k]
   end

   # Calculate divergence of U
   U = state_GV[:,:,:,U_]./state_GV[:,:,:,Rho_] # This costs memory, but ...
   divergence!(param, U, div_G)

   if !param.UseConservative
      @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI
         source_GV[i,j,k,P_] = -(γ-1.0)*state_GV[i+nG,j+nG,k+nG,P_]*div_G[i,j,k]
      end
   else
      @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
         ub = state_GV[i,j,k,Ux_]*state_GV[i,j,k,Bx_] +
            state_GV[i,j,k,Uy_]*state_GV[i,j,k,By_] +
            state_GV[i,j,k,Uz_]*state_GV[i,j,k,Bz_]
         source_GV[i-nG,j-nG,k-nG,E_] = -ub*div_G[i-nG,j-nG,k-nG]
      end
   end

   return source_GV
end

end
