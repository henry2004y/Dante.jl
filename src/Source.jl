module Source

export calc_source!, init_source

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Mx_, My_, Mz_, Bx_, By_, Bz_
using ..Parameters: P_, E_, U_, M_, B_, γ
using ..Divergence: divergence!

function init_source(param::Param)
   GridSize, FullSize = param.GridSize, param.FullSize
   nVar = param.nVar
   source_GV = Array{Float64,4}(undef, GridSize..., nVar)
   div_G = Array{Float64,3}(undef, GridSize...)

   U = Array{Float64,4}(undef, FullSize..., 3)

   return source_GV, U, div_G
end

function calc_source!(param::Param, state_GV, source_GV, U, div_G)

   nI, nJ, nK, nG = param.nI, param.nJ, param.nK, param.nG
   iMin, iMax, jMin, jMax, kMin, kMax =
      param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   # Compute ∇⋅B using central difference
   @views divergence!(param, state_GV[:,:,:,B_], div_G)

   source_GV[:,:,:,Rho_] .= 0.0

   @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI
      source_GV[i,j,k,Mx_] = -state_GV[i+nG,j+nG,k+nG,Bx_]*div_G[i,j,k]
      source_GV[i,j,k,My_] = -state_GV[i+nG,j+nG,k+nG,By_]*div_G[i,j,k]
      source_GV[i,j,k,Mz_] = -state_GV[i+nG,j+nG,k+nG,Bz_]*div_G[i,j,k]
      source_GV[i,j,k,Bx_] =  state_GV[i+nG,j+nG,k+nG,Ux_]*div_G[i,j,k]
      source_GV[i,j,k,By_] =  state_GV[i+nG,j+nG,k+nG,Uy_]*div_G[i,j,k]
      source_GV[i,j,k,Bz_] =  state_GV[i+nG,j+nG,k+nG,Uz_]*div_G[i,j,k]
   end

   # Calculate divergence of U
   for k = 1:size(U,3), j = 1:size(U,2), i = 1:size(U,1)
      U[i,j,k,1] = state_GV[i,j,k,Mx_] / state_GV[i,j,k,Rho_]
      U[i,j,k,2] = state_GV[i,j,k,My_] / state_GV[i,j,k,Rho_]
      U[i,j,k,3] = state_GV[i,j,k,Mz_] / state_GV[i,j,k,Rho_]
   end
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
