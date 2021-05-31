module FaceValue

export calc_face_value!, init_face_value, FaceState, FaceGradient

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_

abstract type FaceState end

struct FaceStateCopy{T<:AbstractFloat} <: FaceState
   LState_XV::Array{T,4}
   RState_XV::Array{T,4}
   LState_YV::Array{T,4}
   RState_YV::Array{T,4}
   LState_ZV::Array{T,4}
   RState_ZV::Array{T,4}
end

struct FaceStateView{T<:AbstractFloat} <: FaceState
   LState_XV::SubArray{T,4,Array{T,4},
      Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64},
      Base.Slice{Base.OneTo{Int64}}},false}
   RState_XV::SubArray{T,4,Array{T,4},
      Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64},
      Base.Slice{Base.OneTo{Int64}}},false}
   LState_YV::SubArray{T,4,Array{T,4},
      Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64},
      Base.Slice{Base.OneTo{Int64}}},false}
   RState_YV::SubArray{T,4,Array{T,4},
      Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64},
      Base.Slice{Base.OneTo{Int64}}},false}
   LState_ZV::SubArray{T,4,Array{T,4},
      Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64},
      Base.Slice{Base.OneTo{Int64}}},false}
   RState_ZV::SubArray{T,4,Array{T,4},
      Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64},
      Base.Slice{Base.OneTo{Int64}}},false}
end

struct FaceGradient{T<:AbstractFloat}
   dq_X::Array{T,4}
   dq_Y::Array{T,4}
   dq_Z::Array{T,4}
end

function init_face_value(param::Param)

   nI, nJ, nK, nVar = param.nI, param.nJ, param.nK, param.nVar

   if param.Order == 2
      dq_X = Array{Float64,4}(undef,nI+2,nJ,nK,nVar)
      dq_Y = Array{Float64,4}(undef,nI,nJ+2,nK,nVar)
      dq_Z = Array{Float64,4}(undef,nI,nJ,nK+2,nVar)

      faceGradient = FaceGradient(dq_X, dq_Y, dq_Z)

      LState_XV = Array{Float64,4}(undef,1+nI,nJ,nK,nVar)
      RState_XV = Array{Float64,4}(undef,1+nI,nJ,nK,nVar)
      LState_YV = Array{Float64,4}(undef,nI,nJ+1,nK,nVar)
      RState_YV = Array{Float64,4}(undef,nI,nJ+1,nK,nVar)
      LState_ZV = Array{Float64,4}(undef,nI,nJ,nK+1,nVar)
      RState_ZV = Array{Float64,4}(undef,nI,nJ,nK+1,nVar)

      faceState = FaceStateCopy(LState_XV,RState_XV,LState_YV,RState_YV,
         LState_ZV,RState_ZV)
   elseif param.Order == 1
      # Fake pointers
      dq_X = Array{Float64,4}(undef,0,0,0,0)
      dq_Y = Array{Float64,4}(undef,0,0,0,0)
      dq_Z = Array{Float64,4}(undef,0,0,0,0)

      temp = Array{Float64,4}(undef,1,1,1,1)
      LState_XV = @view temp[1:1,1:1,1:1,:]
      RState_XV = @view temp[1:1,1:1,1:1,:]
      LState_YV = @view temp[1:1,1:1,1:1,:]
      RState_YV = @view temp[1:1,1:1,1:1,:]
      LState_ZV = @view temp[1:1,1:1,1:1,:]
      RState_ZV = @view temp[1:1,1:1,1:1,:]

      faceGradient = FaceGradient(dq_X, dq_Y, dq_Z)
      faceState = FaceStateView(LState_XV, RState_XV, LState_YV, RState_YV,
         LState_ZV, RState_ZV)
   end

   return faceState, faceGradient
end

"""
	calc_face_value!(param, state_GV, faceState, faceGradient)

Type instability for the return type is introduced, but it seems ok.
I don't know how to modify faceState for views in 1st order.
"""
function calc_face_value!(param::Param, state_GV,
   faceState::FaceState, faceGradient::FaceGradient)

   if param.Order == 1
      iMin, iMax, jMin, jMax, kMin, kMax =
      param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

      LState_XV = @view state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
      RState_XV = @view state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:]
      LState_YV = @view state_GV[iMin:iMax,jMin-1:jMax,kMin:kMax,:]
      RState_YV = @view state_GV[iMin:iMax,jMin:jMax+1,kMin:kMax,:]
      LState_ZV = @view state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax,:]
      RState_ZV = @view state_GV[iMin:iMax,jMin:jMax,kMin:kMax+1,:]

      faceState = FaceStateView(LState_XV, RState_XV, LState_YV, RState_YV,
         LState_ZV, RState_ZV)
   elseif param.Order == 2
      # Compute and limit slopes
      nI,nJ,nK,nG,nVar = param.nI, param.nJ, param.nK, param.nG, param.nVar

      LState_XV, RState_XV = faceState.LState_XV, faceState.RState_XV
      LState_YV, RState_YV = faceState.LState_YV, faceState.RState_YV
      LState_ZV, RState_ZV = faceState.LState_ZV, faceState.RState_ZV

      dq_X = faceGradient.dq_X
      dq_Y = faceGradient.dq_Y
      dq_Z = faceGradient.dq_Z

      # Get slope with limiters
      if param.limiter == "MC"
         # Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}

         for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI+2
            dqR_X = state_GV[i+2,j+nG,k+nG,iVar] - state_GV[i+1,j+nG,k+nG,iVar]
            dqL_X = state_GV[i+1,j+nG,k+nG,iVar] - state_GV[i,j+nG,k+nG,iVar]
            dqC_X = state_GV[i+2,j+nG,k+nG,iVar] - state_GV[i,j+nG,k+nG,iVar]
            dq_X[i,j,k,iVar] = minmod(dqR_X, dqL_X, dqC_X)
         end

         for iVar=1:nVar, k=1:nK, j=1:nJ+2, i=1:nI
            dqR_Y = state_GV[i+nG,j+2,k+nG,iVar] - state_GV[i+nG,j+1,k+nG,iVar]
            dqL_Y = state_GV[i+nG,j+1,k+nG,iVar] - state_GV[i+nG,j,k+nG,iVar]
            dqC_Y = state_GV[i+nG,j+2,k+nG,iVar] - state_GV[i+nG,j,k+nG,iVar]
            dq_Y[i,j,k,iVar] = minmod(dqR_Y, dqL_Y, dqC_Y)
         end

         for iVar=1:nVar, k=1:nK+2, j=1:nJ, i=1:nI
            dqR_Z = state_GV[i+nG,j+nG,k+2,iVar] - state_GV[i+nG,j+nG,k+1,iVar]
            dqL_Z = state_GV[i+nG,j+nG,k+1,iVar] - state_GV[i+nG,j+nG,k,iVar]
            dqC_Y = state_GV[i+nG,j+nG,k+2,iVar] - state_GV[i+nG,j+nG,k,iVar]
            dq_Z[i,j,k,iVar] = minmod(dqR_Z, dqL_Z, dqL_Z)
         end

      elseif param.limiter == "MM" # Minmod limiter
         # Find dq_j = minmod{fwd diff, bwd diff}

         for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI+2
            dqR_X = state_GV[i+2,j+nG,k+nG,iVar] - state_GV[i+1,j+nG,k+nG,iVar]
            dqL_X = state_GV[i+1,j+nG,k+nG,iVar] - state_GV[i,j+nG,k+nG,iVar]
            dq_X[i,j,k,iVar] = minmod(dqR_X, dqL_X)
         end

         for iVar=1:nVar, k=1:nK, j=1:nJ+2, i=1:nI
            dqR_Y = state_GV[i+nG,j+2,k+nG,iVar] - state_GV[i+nG,j+1,k+nG,iVar]
            dqL_Y = state_GV[i+nG,j+1,k+nG,iVar] - state_GV[i+nG,j,k+nG,iVar]
            dq_Y[i,j,k,iVar] = minmod(dqR_Y, dqL_Y)
         end

         for iVar=1:nVar, k=1:nK+2, j=1:nJ, i=1:nI
            dqR_Z = state_GV[i+nG,j+nG,k+2,iVar] - state_GV[i+nG,j+nG,k+1,iVar]
            dqL_Z = state_GV[i+nG,j+nG,k+1,iVar] - state_GV[i+nG,j+nG,k,iVar]
            dq_Z[i,j,k,iVar] = minmod(dqR_Z, dqL_Z)
         end
      end

      # Linear interpolation onto edge centers
      for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
         LState_XV[i,j,k,iVar] = state_GV[i+nG-1,j+nG,k+nG,iVar] +
            0.5*dq_X[i,j,k,iVar]
      end
      for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
         RState_XV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG,iVar] -
            0.5*dq_X[i+1,j,k,iVar]
      end
      for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
         LState_YV[i,j,k,iVar] = state_GV[i+nG,j+nG-1,k+nG,iVar] +
            0.5*dq_Y[i,j,k,iVar]
      end
      for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
         RState_YV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG,iVar] -
            0.5*dq_Y[i,j+1,k,iVar]
      end
      for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
         LState_ZV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG-1,iVar] +
            0.5*dq_Z[i,j,k,iVar]
      end
      for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
         RState_ZV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG,iVar] -
            0.5*dq_Z[i,j,k+1,iVar]
      end
   end

   return faceState
end

"""
	minmod(a, b)
OUTPUT:
 m: zero if opposite sign, otherwise the one of smaller magnitude.
"""
minmod(a, b) = (sign(a) + sign(b)) / 2.0 * min(abs(a), abs(b))

"""
	minmod(a, b, c)
For three inputs, use Harten's generalized definition.
OUTPUT:
 m: zero if opposite sign, otherwise the one of smaller magnitude.
"""
function minmod(a, b, c)
   s = (sign(a) + sign(b) + sign(c))/3.0
   if abs(s) == 1
      m = s*min(abs(a),abs(b),abs(c))
   else
      m = 0.0
   end
   return m
end

end
