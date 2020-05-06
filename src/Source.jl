module Source

export calc_source!, init_source, Div

using LoopVectorization

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Mx_, My_, Mz_, Bx_, By_, Bz_
using ..Parameters: P_, E_, U_, M_, B_, γ
using ..Divergence: divergence!

struct Div{T<:AbstractFloat}
   div_G::Array{T,3}
end

function init_source(param::Param)
   GridSize, FullSize = param.GridSize, param.FullSize
   nVar = param.nVar
   source_GV = Array{Float64,4}(undef, GridSize..., nVar)
   div_G = Array{Float64,3}(undef, GridSize...)

   div = Div(div_G)

   U = Array{Float64,4}(undef, FullSize..., 3)

   return source_GV, U, div
end

function calc_source!(param::Param, state_GV, source_GV, U::Array{Float64,4},
   div::Div)

   nVar, nG = param.nVar, param.nG
   x, y, z = param.x, param.y, param.z

   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   # I can preallocate div_G and DivU!
   nI,nJ,nK = param.nI, param.nJ, param.nK

   div_G = div.div_G
   # Calculate divergence of B using central difference
   B = @view state_GV[:,:,:,B_]
   divergence!(param, B, div_G)

   source_GV[:,:,:,Rho_] .= 0.0

   @avx for k = 1:nK, j = 1:nJ, i = 1:nI
      source_GV[i,j,k,Mx_] = -state_GV[i+nG,j+nG,k+nG,Bx_]*div_G[i,j,k]
      source_GV[i,j,k,My_] = -state_GV[i+nG,j+nG,k+nG,By_]*div_G[i,j,k]
      source_GV[i,j,k,Mz_] = -state_GV[i+nG,j+nG,k+nG,Bz_]*div_G[i,j,k]
      source_GV[i,j,k,Bx_] =  state_GV[i+nG,j+nG,k+nG,Ux_]*div_G[i,j,k]
      source_GV[i,j,k,By_] =  state_GV[i+nG,j+nG,k+nG,Uy_]*div_G[i,j,k]
      source_GV[i,j,k,Bz_] =  state_GV[i+nG,j+nG,k+nG,Uz_]*div_G[i,j,k]
   end

   # Calculate divergence of U
   @. U = state_GV[:,:,:,M_] / state_GV[:,:,:,Rho_]
   divergence!(param, U, div_G)

   if !param.UseConservative
      @avx for k = 1:nK, j = 1:nJ, i = 1:nI
         source_GV[i,j,k,P_] = -(γ-1.0)*state_GV[i+nG,j+nG,k+nG,P_]*div_G[i,j,k]
      end
   else
      @avx for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
         ub = state_GV[i,j,k,Ux_]*state_GV[i,j,k,Bx_] +
            state_GV[i,j,k,Uy_]*state_GV[i,j,k,By_] +
            state_GV[i,j,k,Uz_]*state_GV[i,j,k,Bz_]
         source_GV[i-nG,j-nG,k-nG,E_] = -ub*div_G[i-nG,j-nG,k-nG]
      end
   end

   return source_GV
end

end
