# Face value

struct FaceState{T<:AbstractArray}
   LState_XV::T
   RState_XV::T
   LState_YV::T
   RState_YV::T
   LState_ZV::T
   RState_ZV::T
end

struct FaceGradient{T<:AbstractFloat}
   dq_X::Array{T,4}
   dq_Y::Array{T,4}
   dq_Z::Array{T,4}
end

function init_face_value(param::Param)

   @unpack nI, nJ, nK, nVar = param

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

      faceState = FaceState(LState_XV,RState_XV,LState_YV,RState_YV,LState_ZV,RState_ZV)
   elseif param.Order == 1
      # Fake pointers
      dq_X = Array{Float64,4}(undef,0,0,0,0)
      dq_Y = Array{Float64,4}(undef,0,0,0,0)
      dq_Z = Array{Float64,4}(undef,0,0,0,0)

      FullSize = param.FullSize
      nVar = param.nVar
   
      temp = zeros(FullSize[1],FullSize[2],FullSize[3], nVar)

      iMin, iMax, jMin, jMax, kMin, kMax =
         param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

      LState_XV = @view temp[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
      RState_XV = @view temp[iMin:iMax+1,jMin:jMax,kMin:kMax,:]
      LState_YV = @view temp[iMin:iMax,jMin-1:jMax,kMin:kMax,:]
      RState_YV = @view temp[iMin:iMax,jMin:jMax+1,kMin:kMax,:]
      LState_ZV = @view temp[iMin:iMax,jMin:jMax,kMin-1:kMax,:]
      RState_ZV = @view temp[iMin:iMax,jMin:jMax,kMin:kMax+1,:]

      faceGradient = FaceGradient(dq_X, dq_Y, dq_Z)
      faceState = FaceState(LState_XV, RState_XV, LState_YV, RState_YV,LState_ZV, RState_ZV)
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
      @unpack iMin, iMax, jMin, jMax, kMin, kMax = param

      faceState.LState_XV .= @view state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
      faceState.RState_XV .= @view state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:]
      faceState.LState_YV .= @view state_GV[iMin:iMax,jMin-1:jMax,kMin:kMax,:]
      faceState.RState_YV .= @view state_GV[iMin:iMax,jMin:jMax+1,kMin:kMax,:]
      faceState.LState_ZV .= @view state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax,:]
      faceState.RState_ZV .= @view state_GV[iMin:iMax,jMin:jMax,kMin:kMax+1,:]

   elseif param.Order == 2
      # Compute and limit slopes
      @unpack nI, nJ, nK, nG, nVar = param

      LState_XV, RState_XV = faceState.LState_XV, faceState.RState_XV
      LState_YV, RState_YV = faceState.LState_YV, faceState.RState_YV
      LState_ZV, RState_ZV = faceState.LState_ZV, faceState.RState_ZV

      dq_X = faceGradient.dq_X
      dq_Y = faceGradient.dq_Y
      dq_Z = faceGradient.dq_Z

      # Get slope with limiters
      if param.limiter == "MC"
         # Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}

         @inbounds for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI+2
            dqR_X = state_GV[i+2,j+nG,k+nG,iVar] - state_GV[i+1,j+nG,k+nG,iVar]
            dqL_X = state_GV[i+1,j+nG,k+nG,iVar] - state_GV[i,j+nG,k+nG,iVar]
            dqC_X = state_GV[i+2,j+nG,k+nG,iVar] - state_GV[i,j+nG,k+nG,iVar]
            dq_X[i,j,k,iVar] = minmod(dqR_X, dqL_X, dqC_X)
         end

         @inbounds for iVar=1:nVar, k=1:nK, j=1:nJ+2, i=1:nI
            dqR_Y = state_GV[i+nG,j+2,k+nG,iVar] - state_GV[i+nG,j+1,k+nG,iVar]
            dqL_Y = state_GV[i+nG,j+1,k+nG,iVar] - state_GV[i+nG,j,k+nG,iVar]
            dqC_Y = state_GV[i+nG,j+2,k+nG,iVar] - state_GV[i+nG,j,k+nG,iVar]
            dq_Y[i,j,k,iVar] = minmod(dqR_Y, dqL_Y, dqC_Y)
         end

         @inbounds for iVar=1:nVar, k=1:nK+2, j=1:nJ, i=1:nI
            dqR_Z = state_GV[i+nG,j+nG,k+2,iVar] - state_GV[i+nG,j+nG,k+1,iVar]
            dqL_Z = state_GV[i+nG,j+nG,k+1,iVar] - state_GV[i+nG,j+nG,k,iVar]
            dqC_Z = state_GV[i+nG,j+nG,k+2,iVar] - state_GV[i+nG,j+nG,k,iVar]
            dq_Z[i,j,k,iVar] = minmod(dqR_Z, dqL_Z, dqC_Z)
         end

      elseif param.limiter == "MM" # Minmod limiter
         # Find dq_j = minmod{fwd diff, bwd diff}

         @inbounds for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI+2
            dqR_X = state_GV[i+2,j+nG,k+nG,iVar] - state_GV[i+1,j+nG,k+nG,iVar]
            dqL_X = state_GV[i+1,j+nG,k+nG,iVar] - state_GV[i,j+nG,k+nG,iVar]
            dq_X[i,j,k,iVar] = minmod(dqR_X, dqL_X)
         end

         @inbounds for iVar=1:nVar, k=1:nK, j=1:nJ+2, i=1:nI
            dqR_Y = state_GV[i+nG,j+2,k+nG,iVar] - state_GV[i+nG,j+1,k+nG,iVar]
            dqL_Y = state_GV[i+nG,j+1,k+nG,iVar] - state_GV[i+nG,j,k+nG,iVar]
            dq_Y[i,j,k,iVar] = minmod(dqR_Y, dqL_Y)
         end

         @inbounds for iVar=1:nVar, k=1:nK+2, j=1:nJ, i=1:nI
            dqR_Z = state_GV[i+nG,j+nG,k+2,iVar] - state_GV[i+nG,j+nG,k+1,iVar]
            dqL_Z = state_GV[i+nG,j+nG,k+1,iVar] - state_GV[i+nG,j+nG,k,iVar]
            dq_Z[i,j,k,iVar] = minmod(dqR_Z, dqL_Z)
         end
      end

      # Linear interpolation onto edge centers
      @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
         LState_XV[i,j,k,iVar] = state_GV[i+nG-1,j+nG,k+nG,iVar] + 0.5*dq_X[i,j,k,iVar]
      end
      @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
         RState_XV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG,iVar] - 0.5*dq_X[i+1,j,k,iVar]
      end
      @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
         LState_YV[i,j,k,iVar] = state_GV[i+nG,j+nG-1,k+nG,iVar] + 0.5*dq_Y[i,j,k,iVar]
      end
      @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
         RState_YV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG,iVar] - 0.5*dq_Y[i,j+1,k,iVar]
      end
      @inbounds for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
         LState_ZV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG-1,iVar] + 0.5*dq_Z[i,j,k,iVar]
      end
      @inbounds for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
         RState_ZV[i,j,k,iVar] = state_GV[i+nG,j+nG,k+nG,iVar] - 0.5*dq_Z[i,j,k+1,iVar]
      end
   end
   return
end

"""
	minmod(a, b)

Return zero if opposite sign, otherwise the one of smaller magnitude.
"""
minmod(a, b) = (sign(a) + sign(b)) / 2.0 * min(abs(a), abs(b))

"""
	minmod(a, b, c)

For three inputs, use Harten's generalized definition.
Return zero if opposite sign, otherwise the one of smaller magnitude.
"""
function minmod(a, b, c)
   s = (sign(a) + sign(b) + sign(c)) / 3.0
   if abs(s) == 1
      m = s*min(abs(a), abs(b), abs(c))
   else
      m = 0.0
   end
   m
end