module FaceValue

export calc_face_value, FaceState

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_

struct FaceState
   LState_XV::Array{Float64,4}
   RState_XV::Array{Float64,4}
   LState_YV::Array{Float64,4}
   RState_YV::Array{Float64,4}
   LState_ZV::Array{Float64,4}
   RState_ZV::Array{Float64,4}
end

function calc_face_value(param::Param, state_GV::Array{Float64,4})

   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   if param.Order == 1
      LState_XV = state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
      RState_XV = state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:]
      LState_YV = state_GV[iMin:iMax,jMin-1:jMax,kMin:kMax,:]
      RState_YV = state_GV[iMin:iMax,jMin:jMax+1,kMin:kMax,:]
      LState_ZV = state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax,:]
      RState_ZV = state_GV[iMin:iMax,jMin:jMax,kMin:kMax+1,:]
   elseif param.Order == 2
      # Compute and limit slopes

      # Get slope with limiters
      if param.limiter == "MC"
         # Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
         dqR_X =
            state_GV[iMin:iMax+2,  jMin:jMax,kMin:kMax,:] .-
            state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:]
         dqL_X =
            state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:] .-
            state_GV[iMin-2:iMax  ,jMin:jMax,kMin:kMax,:]
         dqC_X =
            state_GV[iMin:iMax+2,  jMin:jMax,kMin:kMax,:] .-
            state_GV[iMin-2:iMax  ,jMin:jMax,kMin:kMax,:]

         dq_X = minmod.(dqR_X,dqL_X,dqC_X)

         dqR_Y =
            state_GV[iMin:iMax,jMin:jMax+2,  kMin:kMax,:] .-
            state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:]
         dqL_Y =
            state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:] .-
            state_GV[iMin:iMax,jMin-2:jMax  ,kMin:kMax,:]

         dq_Y = minmod.(dqR_Y,dqL_Y,dqC_Y)

         dqR_Z =
            state_GV[iMin:iMax,jMin:jMax,kMin:kMax+2,  :] .-
            state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:]
         dqL_Z =
            state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:] .-
            state_GV[iMin:iMax,jMin:jMax,kMin-2:kMax  ,:]

         dq_Z = minmod.(dqR_Z,dqL_Z,dqC_Z)
      elseif param.limiter == "MM" # Minmod limiter
         # Find dq_j = minmod{fwd diff, bwd diff}

         dqR_X =
            state_GV[iMin:iMax+2  ,jMin:jMax,kMin:kMax,:] .-
            state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:]
         dqL_X =
            state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:] .-
            state_GV[iMin-2:iMax  ,jMin:jMax,kMin:kMax,:]

         dq_X = minmod.(dqR_X,dqL_X)

         dqR_Y =
            state_GV[iMin:iMax,jMin:jMax+2,  kMin:kMax,:] .-
            state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:]
         dqL_Y =
            state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:] .-
            state_GV[iMin:iMax,jMin-2:jMax  ,kMin:kMax,:]

         dq_Y = minmod.(dqR_Y,dqL_Y)

         dqR_Z =
            state_GV[iMin:iMax,jMin:jMax,kMin:kMax+2  ,:] .-
            state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:]
         dqL_Z =
            state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:] .-
            state_GV[iMin:iMax,jMin:jMax,kMin-2:kMax  ,:]

         dq_Z = minmod.(dqR_Z,dqL_Z)
      end

      # Linear interpolation onto edge centers
      LState_XV = state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:] .+
         0.5*dq_X[1:end-1,:,:,:]
      RState_XV = state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:] .-
         0.5*dq_X[2:end  ,:,:,:]
      LState_YV = state_GV[iMin:iMax,jMin-1:jMax,kMin:kMax,:] .+
         0.5*dq_Y[:,1:end-1,:,:]
      RState_YV = state_GV[iMin:iMax,jMin:jMax+1,kMin:kMax,:] .-
         0.5*dq_Y[:,2:end  ,:,:]
      LState_ZV = state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax,:] .+
         0.5*dq_Z[:,:,1:end-1,:]
      RState_ZV = state_GV[iMin:iMax,jMin:jMax,kMin:kMax+1,:] .-
         0.5*dq_Z[:,:,2:end  ,:]

   end

   faceValue = FaceState(LState_XV,RState_XV,LState_YV,RState_YV,LState_ZV,RState_ZV)

   return faceValue
end

"""
   minmod(a,b)
OUTPUT:
 m: zero if opposite sign, otherwise the one of smaller magnitude.
"""
minmod(a::Float64,b::Float64) = (sign(a) + sign(b)) / 2.0 * min(abs(a), abs(b))
#function minmod(a,b)
#   m = @. (sign(a) + sign(b)) / 2.0 * min(abs(a), abs(b))
#end

"""
   minmod(a,b,c)
For three inputs, use Harten's generalized definition.
OUTPUT:
 m: zero if opposite sign, otherwise the one of smaller magnitude.
"""
function minmod(a::Float64,b::Float64,c::Float64)
   s = (sign(a) + sign(b) + sign(c))/3.0
   if abs(s) == 1
      m = s*min(abs(a),abs(b),abs(c))
   else
      m = 0.0
   end
   return m
end

end
