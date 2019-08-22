module Source

export calc_source

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_
using ..Parameters: γ
using ..Divergence: divergence_ndgrid

function calc_source(param::Param, state_GV::Array{Float64,4})

   GridSize = param.GridSize
   nI, nJ, nK = param.nI, param.nJ, param.nK
   nVar, nG = param.nVar, param.nG
   x, y, z = param.x, param.y, param.z

   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   source_GV = Array{Float64,4}(undef, GridSize..., nVar)

   # Calculate divergence of B using central difference
   DivB = divergence_ndgrid(x,y,z,state_GV[:,:,:,B_]) # takes time

   # Calculate divergence of U
   DivU = divergence_ndgrid(x,y,z,state_GV[:,:,:,U_]./state_GV[:,:,:,Rho_]) # takes time

   source_GV[:,:,:,Rho_] .= 0.0

   @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
      source_GV[i-nG,j-nG,k-nG,Ux_] = -state_GV[i,j,k,Bx_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Uy_] = -state_GV[i,j,k,By_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Uz_] = -state_GV[i,j,k,Bz_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Bx_] =  state_GV[i,j,k,Ux_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,By_] =  state_GV[i,j,k,Uy_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Bz_] =  state_GV[i,j,k,Uz_]*DivB[i,j,k]
   end

   if !param.UseConservative
      @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
         source_GV[i-nG,j-nG,k-nG,P_] = -(γ-1.0)*state_GV[i,j,k,P_]*DivU[i,j,k]
      end
   else
      @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
         ub = state_GV[i,j,k,Ux_]*state_GV[i,j,k,Bx_] +
            state_GV[i,j,k,Uy_]*state_GV[i,j,k,By_] +
            state_GV[i,j,k,Uz_]*state_GV[i,j,k,Bz_]
         source_GV[i-nG,j-nG,k-nG,E_] = -ub*DivB[i,j,k]
      end
   end

   return source_GV
end

end
