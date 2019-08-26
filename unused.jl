# unused

# I feel that dot fusion may be the culprit of slowness.
# Let me try loops instead.

#iMin,iMax,jMin,jMax,kMin,kMax = param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

# face value
LState_XV::SubArray
RState_XV::SubArray
LState_YV::SubArray
RState_YV::SubArray
LState_ZV::SubArray
RState_ZV::SubArray

LState_XV = @view state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
RState_XV = @view state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:]
LState_YV = @view state_GV[iMin:iMax,jMin-1:jMax,kMin:kMax,:]
RState_YV = @view state_GV[iMin:iMax,jMin:jMax+1,kMin:kMax,:]
LState_ZV = @view state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax,:]
RState_ZV = @view state_GV[iMin:iMax,jMin:jMax,kMin:kMax+1,:]

# Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
@views dqR_X =
   state_GV[iMin:iMax+2,  jMin:jMax,kMin:kMax,:] .-
   state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:]
@views dqL_X =
   state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:] .-
   state_GV[iMin-2:iMax  ,jMin:jMax,kMin:kMax,:]
@views dqC_X =
   state_GV[iMin:iMax+2,  jMin:jMax,kMin:kMax,:] .-
   state_GV[iMin-2:iMax  ,jMin:jMax,kMin:kMax,:]

dq_X = minmod.(dqR_X,dqL_X,dqC_X)

@views dqR_Y =
   state_GV[iMin:iMax,jMin:jMax+2,  kMin:kMax,:] .-
   state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:]
@views dqL_Y =
   state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:] .-
   state_GV[iMin:iMax,jMin-2:jMax  ,kMin:kMax,:]
@views dqC_Y =
   state_GV[iMin:iMax,jMin:jMax+2  ,kMin:kMax,:] .-
   state_GV[iMin:iMax,jMin-2:jMax  ,kMin:kMax,:]

dq_Y = minmod.(dqR_Y,dqL_Y,dqC_Y)

@views dqR_Z =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax+2,  :] .-
   state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:]
@views dqL_Z =
   state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:] .-
   state_GV[iMin:iMax,jMin:jMax,kMin-2:kMax  ,:]
@views dqC_Z =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax+2,  :] .-
   state_GV[iMin:iMax,jMin:jMax,kMin-2:kMax,  :]

dq_Z = minmod.(dqR_Z,dqL_Z,dqC_Z)

# Find dq_j = minmod{fwd diff, bwd diff}

@views dqR_X =
   state_GV[iMin:iMax+2  ,jMin:jMax,kMin:kMax,:] .-
   state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:]
@views dqL_X =
   state_GV[iMin-1:iMax+1,jMin:jMax,kMin:kMax,:] .-
   state_GV[iMin-2:iMax  ,jMin:jMax,kMin:kMax,:]

dq_X = minmod.(dqR_X, dqL_X)

@views dqR_Y =
   state_GV[iMin:iMax,jMin:jMax+2,  kMin:kMax,:] .-
   state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:]
@views dqL_Y =
   state_GV[iMin:iMax,jMin-1:jMax+1,kMin:kMax,:] .-
   state_GV[iMin:iMax,jMin-2:jMax  ,kMin:kMax,:]

dq_Y = minmod.(dqR_Y, dqL_Y)

@views dqR_Z =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax+2  ,:] .-
   state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:]
@views dqL_Z =
   state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax+1,:] .-
   state_GV[iMin:iMax,jMin:jMax,kMin-2:kMax  ,:]

dq_Z = minmod.(dqR_Z, dqL_Z)

# Linear interpolation onto edge centers
@views LState_XV = state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:] .+
   0.5*dq_X[1:end-1,:,:,:]
@views RState_XV = state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:] .-
   0.5*dq_X[2:end  ,:,:,:]
@views LState_YV = state_GV[iMin:iMax,jMin-1:jMax,kMin:kMax,:] .+
   0.5*dq_Y[:,1:end-1,:,:]
@views RState_YV = state_GV[iMin:iMax,jMin:jMax+1,kMin:kMax,:] .-
   0.5*dq_Y[:,2:end  ,:,:]
@views LState_ZV = state_GV[iMin:iMax,jMin:jMax,kMin-1:kMax,:] .+
   0.5*dq_Z[:,:,1:end-1,:]
@views RState_ZV = state_GV[iMin:iMax,jMin:jMax,kMin:kMax+1,:] .-
   0.5*dq_Z[:,:,2:end  ,:]


# face flux
bL = zeros(nI+1,nJ+1,nK+1)
bR = zeros(nI+1,nJ+1,nK+1)

bL = zeros(nI+1,nJ+1,nK+1)
bR = zeros(nI+1,nJ+1,nK+1)
uL = zeros(nI+1,nJ+1,nK+1)
uR = zeros(nI+1,nJ+1,nK+1)

bL = sum(LState_XV[:,:,:,B_].^2,dims=4)
bR = sum(RState_XV[:,:,:,B_].^2,dims=4)

bL = sum(LState_YV[:,:,:,B_].^2,dims=4) # takes time
bR = sum(RState_YV[:,:,:,B_].^2,dims=4)

bL = sum(LState_ZV[:,:,:,B_].^2,dims=4)
bR = sum(RState_ZV[:,:,:,B_].^2,dims=4) # takes time

@. LFlux_XV[:,:,:,Rho_] = LState_XV[:,:,:,Ux_]
@. RFlux_XV[:,:,:,Rho_] = RState_XV[:,:,:,Ux_]
@. LFlux_YV[:,:,:,Rho_] = LState_YV[:,:,:,Uy_]
@. RFlux_YV[:,:,:,Rho_] = RState_YV[:,:,:,Uy_]
@. LFlux_ZV[:,:,:,Rho_] = LState_ZV[:,:,:,Uz_]
@. RFlux_ZV[:,:,:,Rho_] = RState_ZV[:,:,:,Uz_]

@. LFlux_XV[:,:,:,Ux_] = LState_XV[:,:,:,Ux_]^2 / LState_XV[:,:,:,Rho_] + LState_XV[:,:,:,P_] + 0.5*$dropdims($sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4) - LState_XV[:,:,:,Bx_]^2
@. RFlux_XV[:,:,:,Ux_] = RState_XV[:,:,:,Ux_]^2 / RState_XV[:,:,:,Rho_] + RState_XV[:,:,:,P_] + 0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4) - RState_XV[:,:,:,Bx_]^2
@. LFlux_XV[:,:,:,Uy_] = LState_XV[:,:,:,Ux_] * LState_XV[:,:,:,Uy_] / LState_XV[:,:,:,Rho_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,By_]
@. RFlux_XV[:,:,:,Uy_] = RState_XV[:,:,:,Ux_] * RState_XV[:,:,:,Uy_] / RState_XV[:,:,:,Rho_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,By_]
@. LFlux_XV[:,:,:,Uz_] = LState_XV[:,:,:,Ux_] * LState_XV[:,:,:,Uz_] / LState_XV[:,:,:,Rho_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,Bz_]
@. RFlux_XV[:,:,:,Uz_] = RState_XV[:,:,:,Ux_] * RState_XV[:,:,:,Uz_] / RState_XV[:,:,:,Rho_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,Bz_]

@. LFlux_YV[:,:,:,Ux_] = LState_YV[:,:,:,Uy_] * LState_YV[:,:,:,Ux_] / LState_YV[:,:,:,Rho_] - LState_YV[:,:,:,By_]*LState_YV[:,:,:,Bx_]
@. RFlux_YV[:,:,:,Ux_] = RState_YV[:,:,:,Uy_] * RState_YV[:,:,:,Ux_] / RState_YV[:,:,:,Rho_] - RState_YV[:,:,:,By_]*RState_YV[:,:,:,Bx_]
@. LFlux_YV[:,:,:,Uy_] = LState_YV[:,:,:,Uy_]^2 / LState_YV[:,:,:,Rho_] + LState_YV[:,:,:,P_] + 0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4) - LState_YV[:,:,:,By_]^2
@. RFlux_YV[:,:,:,Uy_] = RState_YV[:,:,:,Uy_]^2 / RState_YV[:,:,:,Rho_] + RState_YV[:,:,:,P_] + 0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4) - RState_YV[:,:,:,By_]^2
@. LFlux_YV[:,:,:,Uz_] = LState_YV[:,:,:,Uy_] * LState_YV[:,:,:,Uz_] / LState_YV[:,:,:,Rho_] - LState_YV[:,:,:,Bx_]*LState_YV[:,:,:,Bz_]
@. RFlux_YV[:,:,:,Uz_] = RState_YV[:,:,:,Uy_] * RState_YV[:,:,:,Uz_] / RState_YV[:,:,:,Rho_] - RState_YV[:,:,:,Bx_]*RState_YV[:,:,:,Bz_]

@. LFlux_ZV[:,:,:,Ux_] = LState_ZV[:,:,:,Uz_] * LState_ZV[:,:,:,Ux_] / LState_ZV[:,:,:,Rho_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,Bx_]
@. RFlux_ZV[:,:,:,Ux_] = RState_ZV[:,:,:,Uz_] * RState_ZV[:,:,:,Ux_] / RState_ZV[:,:,:,Rho_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,Bx_]
@. LFlux_ZV[:,:,:,Uy_] = LState_ZV[:,:,:,Uz_] * LState_ZV[:,:,:,Uy_] / LState_ZV[:,:,:,Rho_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,By_]
@. RFlux_ZV[:,:,:,Uy_] = RState_ZV[:,:,:,Uz_] * RState_ZV[:,:,:,Uy_] / RState_ZV[:,:,:,Rho_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,By_]
@. LFlux_ZV[:,:,:,Uz_] = LState_ZV[:,:,:,Uz_]^2 / LState_ZV[:,:,:,Rho_] + LState_ZV[:,:,:,P_] + 0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4) - LState_ZV[:,:,:,Bz_]^2
@. RFlux_ZV[:,:,:,Uz_] = RState_ZV[:,:,:,Uz_]^2 / RState_ZV[:,:,:,Rho_] + RState_ZV[:,:,:,P_] + 0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4) - RState_ZV[:,:,:,Bz_]^2

@. LFlux_XV[:,:,:,Bx_] = 0.0
@. RFlux_XV[:,:,:,Bx_] = 0.0
@. LFlux_XV[:,:,:,By_] = LState_XV[:,:,:,Ux_]*LState_XV[:,:,:,By_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,Uy_]
@. RFlux_XV[:,:,:,By_] = RState_XV[:,:,:,Ux_]*RState_XV[:,:,:,By_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,Uy_]
@. LFlux_XV[:,:,:,Bz_] = LState_XV[:,:,:,Ux_]*LState_XV[:,:,:,Bz_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,Uz_]
@. RFlux_XV[:,:,:,Bz_] = RState_XV[:,:,:,Ux_]*RState_XV[:,:,:,Bz_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,Uz_]

@. LFlux_YV[:,:,:,Bx_] = LState_YV[:,:,:,Uy_]*LState_YV[:,:,:,Bx_] - LState_YV[:,:,:,By_]*LState_YV[:,:,:,Ux_]
@. RFlux_YV[:,:,:,Bx_] = RState_YV[:,:,:,Uy_]*RState_YV[:,:,:,Bx_] - RState_YV[:,:,:,By_]*RState_YV[:,:,:,Ux_]
@. LFlux_YV[:,:,:,By_] = 0.0
@. RFlux_YV[:,:,:,By_] = 0.0
@. LFlux_YV[:,:,:,Bz_] = LState_YV[:,:,:,Uy_]*LState_YV[:,:,:,Bz_] - LState_YV[:,:,:,By_]*LState_YV[:,:,:,Uz_]
@. RFlux_YV[:,:,:,Bz_] = RState_YV[:,:,:,Uy_]*RState_YV[:,:,:,Bz_] - RState_YV[:,:,:,By_]*RState_YV[:,:,:,Uz_]

@. LFlux_ZV[:,:,:,Bx_] = LState_ZV[:,:,:,Uz_]*LState_ZV[:,:,:,Bx_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,Ux_]
@. RFlux_ZV[:,:,:,Bx_] = RState_ZV[:,:,:,Uz_]*RState_ZV[:,:,:,Bx_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,Ux_]
@. LFlux_ZV[:,:,:,By_] = LState_ZV[:,:,:,Uz_]*LState_ZV[:,:,:,By_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,Uy_]
@. RFlux_ZV[:,:,:,By_] = RState_ZV[:,:,:,Uz_]*RState_ZV[:,:,:,By_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,Uy_]
@. LFlux_ZV[:,:,:,Bz_] = 0.0
@. RFlux_ZV[:,:,:,Bz_] = 0.0

@. LFlux_XV[:,:,:,P_] = LState_XV[:,:,:,P_]*LState_XV[:,:,:,Ux_]/LState_XV[:,:,:,Rho_]
@. RFlux_XV[:,:,:,P_] = RState_XV[:,:,:,P_]*RState_XV[:,:,:,Ux_]/RState_XV[:,:,:,Rho_]
@. LFlux_YV[:,:,:,P_] = LState_YV[:,:,:,P_]*LState_YV[:,:,:,Uy_]/LState_YV[:,:,:,Rho_]
@. RFlux_YV[:,:,:,P_] = RState_YV[:,:,:,P_]*RState_YV[:,:,:,Uy_]/RState_YV[:,:,:,Rho_]
@. LFlux_ZV[:,:,:,P_] = LState_ZV[:,:,:,P_]*LState_ZV[:,:,:,Uz_]/LState_ZV[:,:,:,Rho_]
@. RFlux_ZV[:,:,:,P_] = RState_ZV[:,:,:,P_]*RState_ZV[:,:,:,Uz_]/RState_ZV[:,:,:,Rho_]

@. Flux_XV = 0.5 * (LFlux_XV + RFlux_XV)
@. Flux_YV = 0.5 * (LFlux_YV + RFlux_YV)
@. Flux_ZV = 0.5 * (LFlux_ZV + RFlux_ZV)

u = zeros(nI+1,nJ+1,nK+1)
b = zeros(nI+1,nJ+1,nK+1)
ub= zeros(nI+1,nJ+1,nK+1)

u = sum(LState_XV[:,:,:,U_].^2,dims=4)
b = sum(LState_XV[:,:,:,B_].^2,dims=4)
ub = sum(LState_XV[:,:,:,U_].*LState_XV[:,:,:,B_],dims=4)

@. LFlux_XV[:,:,:,E_] = LState_XV[:,:,:,Ux_]/LState_XV[:,:,:,Rho_]*
      ((LState_XV[:,:,:,P_]/(γ-1.0) + 0.5/LState_XV[:,:,:,Rho_]*
      u + 0.5*b + LState_XV[:,:,:,P_] + 0.5*b)) - ub*LState_XV[:,:,:,Bx_]

u = sum(RState_XV[:,:,:,U_].^2,dims=4) # takes time!
b = sum(RState_XV[:,:,:,B_].^2,dims=4) # takes time!
ub = sum(RState_XV[:,:,:,U_].*RState_XV[:,:,:,B_],dims=4)

@. RFlux_XV[:,:,:,E_] = RState_XV[:,:,:,Ux_]/RState_XV[:,:,:,Rho_]*
      ((RState_XV[:,:,:,P_]/(γ-1.0) + 0.5/RState_XV[:,:,:,Rho_]*
      $dropdims($sum(RState_XV[:,:,:,U_].^2,dims=4);dims=4) +
      0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4) + RState_XV[:,:,:,P_] +
      0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4))) -
      $dropdims($sum(RState_XV[:,:,:,U_].*RState_XV[:,:,:,B_],dims=4);dims=4)*
      RState_XV[:,:,:,Bx_]

u = sum(LState_YV[:,:,:,U_].^2,dims=4)
b = sum(LState_YV[:,:,:,B_].^2,dims=4)
ub = sum(LState_YV[:,:,:,U_].*LState_YV[:,:,:,B_],dims=4)

@. LFlux_YV[:,:,:,E_] = LState_YV[:,:,:,Uy_]/LState_YV[:,:,:,Rho_]*
      ((LState_YV[:,:,:,P_]/(γ-1.0) + 0.5/LState_YV[:,:,:,Rho_]*
      $dropdims($sum(LState_YV[:,:,:,U_].^2,dims=4);dims=4) +
      0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4) + LState_YV[:,:,:,P_] +
      0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4))) -
      $dropdims($sum(LState_YV[:,:,:,U_].*LState_YV[:,:,:,B_],dims=4);dims=4)*
      LState_YV[:,:,:,By_]

u = sum(RState_YV[:,:,:,U_].^2,dims=4) # dropdims takes time
b = sum(RState_YV[:,:,:,B_].^2,dims=4)
ub = sum(RState_YV[:,:,:,U_].*RState_YV[:,:,:,B_],dims=4) # takes time

@. RFlux_YV[:,:,:,E_] = RState_YV[:,:,:,Uy_]/RState_YV[:,:,:,Rho_]*
      ((RState_YV[:,:,:,P_]/(γ-1.0) +
      0.5/RState_YV[:,:,:,Rho_]*$dropdims($sum(RState_YV[:,:,:,U_].^2,dims=4);dims=4) +
      0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4) + RState_YV[:,:,:,P_] +
      0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4))) -
      $dropdims($sum(RState_YV[:,:,:,U_].*RState_YV[:,:,:,B_],dims=4);dims=4)*
      RState_YV[:,:,:,By_]

u = sum(LState_ZV[:,:,:,U_].^2,dims=4)
b = sum(LState_ZV[:,:,:,B_].^2,dims=4)
ub = sum(LState_ZV[:,:,:,U_].*LState_ZV[:,:,:,B_],dims=4)

@. LFlux_ZV[:,:,:,E_] = LState_ZV[:,:,:,Uz_]/LState_ZV[:,:,:,Rho_]*
      ((LState_ZV[:,:,:,P_]/(γ-1.0) + 0.5/LState_ZV[:,:,:,Rho_]*
      $dropdims($sum(LState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
      0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4) + LState_ZV[:,:,:,P_] +
      0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4))) -
      $dropdims($sum(LState_ZV[:,:,:,U_].*LState_ZV[:,:,:,B_],dims=4);dims=4)*
      LState_ZV[:,:,:,Bz_]

u = sum(RState_ZV[:,:,:,U_].^2,dims=4)
b = sum(RState_ZV[:,:,:,B_].^2,dims=4)
ub = sum(RState_ZV[:,:,:,U_].*RState_ZV[:,:,:,B_],dims=4)

@. RFlux_ZV[:,:,:,E_] = RState_ZV[:,:,:,Uz_]/RState_ZV[:,:,:,Rho_]*
   ((RState_ZV[:,:,:,P_]/(γ-1.0) + 0.5/RState_ZV[:,:,:,Rho_]*
   $dropdims($sum(RState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
   0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4) + RState_ZV[:,:,:,P_] +
   0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4))) -
   $dropdims($sum(RState_ZV[:,:,:,U_].*RState_ZV[:,:,:,B_],dims=4);dims=4)*
   RState_ZV[:,:,:,Bz_]



@. Flux_XV -= 0.5*Cmax_XF*(RState_XV - LState_XV)
@. Flux_YV -= 0.5*Cmax_YF*(RState_YV - LState_YV)
@. Flux_ZV -= 0.5*Cmax_ZF*(RState_ZV - LState_ZV)

uL = sum(LState_XV[:,:,:,U_].^2,dims=4)
bL = sum(LState_XV[:,:,:,B_].^2,dims=4)
uR = sum(RState_XV[:,:,:,U_].^2,dims=4)
bR = sum(RState_XV[:,:,:,B_].^2,dims=4)

@. Flux_XV[:,:,:,Rho_:Bz_] -=
0.5*Cmax_XF*(RState_XV[:,:,:,Rho_:Bz_] - LState_XV[:,:,:,Rho_:Bz_])
@. Flux_XV[:,:,:,E_] -= 0.5*Cmax_XF*(
(RState_XV[:,:,:,P_]/(γ-1.0) + 0.5/RState_XV[:,:,:,Rho_]*
$dropdims($sum(RState_XV[:,:,:,U_].^2,dims=4);dims=4) +
0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4)) -
(LState_XV[:,:,:,P_]/(γ-1.0) +
0.5/LState_XV[:,:,:,Rho_]*$dropdims($sum(LState_XV[:,:,:,U_].^2,dims=4);dims=4) +
0.5*$dropdims($sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4)))

uL = sum(LState_YV[:,:,:,U_].^2,dims=4)
bL = sum(LState_YV[:,:,:,B_].^2,dims=4)
uR = sum(RState_YV[:,:,:,U_].^2,dims=4)
bR = sum(RState_YV[:,:,:,B_].^2,dims=4)

@. Flux_YV[:,:,:,Rho_:Bz_] -=
0.5*Cmax_YF*(RState_YV[:,:,:,Rho_:Bz_] - LState_YV[:,:,:,Rho_:Bz_])
@. Flux_YV[:,:,:,E_] -= 0.5*Cmax_YF*(
(RState_YV[:,:,:,P_]/(γ-1.0) + 0.5/RState_YV[:,:,:,Rho_]*
$dropdims($sum(RState_YV[:,:,:,U_].^2,dims=4);dims=4) +
0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4)) -
(LState_YV[:,:,:,P_]/(γ-1.0) +
0.5/LState_YV[:,:,:,Rho_]*$dropdims($sum(LState_YV[:,:,:,U_].^2,dims=4);dims=4) +
0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4)))

uL = sum(LState_ZV[:,:,:,U_].^2,dims=4)
bL = sum(LState_ZV[:,:,:,B_].^2,dims=4)
uR = sum(RState_ZV[:,:,:,U_].^2,dims=4)
bR = sum(RState_ZV[:,:,:,B_].^2,dims=4)

@. Flux_ZV[:,:,:,Rho_:Bz_] -=
0.5*Cmax_ZF*(RState_ZV[:,:,:,Rho_:Bz_] - LState_ZV[:,:,:,Rho_:Bz_])
@. Flux_ZV[:,:,:,E_] -= 0.5*Cmax_ZF*(
(RState_ZV[:,:,:,P_]/(γ-1.0) + 0.5/RState_ZV[:,:,:,Rho_]*
$dropdims($sum(RState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4)) -
(LState_ZV[:,:,:,P_]/(γ-1.0) + 0.5/LState_ZV[:,:,:,Rho_]*
$dropdims($sum(LState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4)))



Cs2_XF = γ.*(LS_XV[:,:,:,P_] .+ RS_XV[:,:,:,P_]) ./
   (LS_XV[:,:,:,Rho_] .+ RS_XV[:,:,:,Rho_])
Cs2_YF = γ.*(LS_YV[:,:,:,P_] .+ RS_YV[:,:,:,P_]) ./
   (LS_YV[:,:,:,Rho_] .+ RS_YV[:,:,:,Rho_])
Cs2_ZF = γ.*(LS_ZV[:,:,:,P_] .+ RS_ZV[:,:,:,P_]) ./
   (LS_ZV[:,:,:,Rho_] .+ RS_ZV[:,:,:,Rho_])

Ca2_XF = ( (LS_XV[:,:,:,Bx_] .+ RS_XV[:,:,:,Bx_]).^2 .+
   (LS_XV[:,:,:,By_] .+ RS_XV[:,:,:,By_]).^2 .+
   (LS_XV[:,:,:,Bz_] .+ RS_XV[:,:,:,Bz_]).^2 ) ./
   (2.0 .*(LS_XV[:,:,:,Rho_] .+ RS_XV[:,:,:,Rho_]))
Ca2_YF = ( (LS_YV[:,:,:,Bx_] .+ RS_YV[:,:,:,Bx_]).^2 .+
   (LS_YV[:,:,:,By_] .+ RS_YV[:,:,:,By_]).^2 .+
   (LS_YV[:,:,:,Bz_] .+ RS_YV[:,:,:,Bz_]).^2 ) ./
   (2.0 .*(LS_YV[:,:,:,Rho_] .+ RS_YV[:,:,:,Rho_]))
Ca2_ZF = ( (LS_ZV[:,:,:,Bx_] .+ RS_ZV[:,:,:,Bx_]).^2 .+
   (LS_ZV[:,:,:,By_] .+ RS_ZV[:,:,:,By_]).^2 .+
   (LS_ZV[:,:,:,Bz_] .+ RS_ZV[:,:,:,Bz_]).^2 ) ./
   (2.0 .*(LS_ZV[:,:,:,Rho_] .+ RS_ZV[:,:,:,Rho_]))

Can2_XF = ( (LS_XV[:,:,:,Bx_] .+ RS_XV[:,:,:,Bx_]).^2 ) ./
   (2.0 .*(LS_XV[:,:,:,Rho_] .+ RS_XV[:,:,:,Rho_]))
Can2_YF = ( (LS_YV[:,:,:,By_] .+ RS_YV[:,:,:,By_]).^2 ) ./
   (2.0 .*(LS_YV[:,:,:,Rho_] + RS_YV[:,:,:,Rho_]))
Can2_ZF = ( (LS_ZV[:,:,:,Bz_] .+ RS_ZV[:,:,:,Bz_]).^2 ) ./
   (2.0 .*(LS_ZV[:,:,:,Rho_] .+ RS_ZV[:,:,:,Rho_]))

@. Cmax_XF = 0.5 * abs(LS_XV[:,:,:,Ux_]/LS_XV[:,:,:,Rho_] +
   RS_XV[:,:,:,Ux_]/RS_XV[:,:,:,Rho_]) + sqrt( 0.5*(Cs2_XF + Ca2_XF +
   sqrt((Cs2_XF + Ca2_XF)^2 - 4.0 *Cs2_XF*Can2_XF)) )

@. Cmax_YF = 0.5 * abs(LS_YV[:,:,:,Uy_]/LS_YV[:,:,:,Rho_] +
   RS_YV[:,:,:,Uy_]/RS_YV[:,:,:,Rho_]) + sqrt( 0.5*(Cs2_YF + Ca2_YF +
   sqrt((Cs2_YF + Ca2_YF)^2 - 4.0 *Cs2_YF*Can2_YF)) )

@. Cmax_ZF = 0.5 * abs(LS_ZV[:,:,:,Uz_]/LS_ZV[:,:,:,Rho_] +
   RS_ZV[:,:,:,Uz_]/RS_ZV[:,:,:,Rho_]) + sqrt( 0.5*(Cs2_ZF + Ca2_ZF +
   sqrt((Cs2_ZF + Ca2_ZF)^2 - 4.0 *Cs2_ZF*Can2_ZF)) )


# get_speed_maxmin
Cs2_LXF = @. γ*LS_XV[:,:,:,P_]/LS_XV[:,:,:,Rho_]
Cs2_LYF = @. γ*LS_YV[:,:,:,P_]/LS_YV[:,:,:,Rho_]
Cs2_LZF = @. γ*LS_ZV[:,:,:,P_]/LS_ZV[:,:,:,Rho_]
Cs2_RXF = @. γ*RS_XV[:,:,:,P_]/RS_XV[:,:,:,Rho_]
Cs2_RYF = @. γ*RS_YV[:,:,:,P_]/RS_YV[:,:,:,Rho_]
Cs2_RZF = @. γ*RS_ZV[:,:,:,P_]/RS_ZV[:,:,:,Rho_]

Ca2_LXF = @. (LS_XV[:,:,:,Bx_]^2 + LS_XV[:,:,:,By_]^2 +
LS_XV[:,:,:,Bz_]^2 ) / LS_XV[:,:,:,Rho_]
Ca2_LYF = @. (LS_YV[:,:,:,Bx_]^2 + LS_YV[:,:,:,By_]^2 +
LS_YV[:,:,:,Bz_]^2 ) / LS_YV[:,:,:,Rho_]
Ca2_LZF = @. (LS_ZV[:,:,:,Bx_]^2 + LS_ZV[:,:,:,By_]^2 +
LS_ZV[:,:,:,Bz_]^2 ) / LS_ZV[:,:,:,Rho_]
Ca2_RXF = @. (RS_XV[:,:,:,Bx_]^2 + RS_XV[:,:,:,By_]^2 +
RS_XV[:,:,:,Bz_]^2 ) / RS_XV[:,:,:,Rho_]
Ca2_RYF = @. (RS_YV[:,:,:,Bx_]^2 + RS_YV[:,:,:,By_]^2 +
RS_YV[:,:,:,Bz_]^2 ) / RS_YV[:,:,:,Rho_]
Ca2_RZF = @. (RS_ZV[:,:,:,Bx_]^2 + RS_ZV[:,:,:,By_]^2 +
RS_ZV[:,:,:,Bz_]^2 ) / RS_ZV[:,:,:,Rho_]


Can2_LXF = @. LS_XV[:,:,:,Bx_]^2 / LS_XV[:,:,:,Rho_]
Can2_LYF = @. LS_YV[:,:,:,By_]^2 / LS_YV[:,:,:,Rho_]
Can2_LZF = @. LS_ZV[:,:,:,Bz_]^2 / LS_ZV[:,:,:,Rho_]
Can2_RXF = @. RS_XV[:,:,:,Bx_]^2 / RS_XV[:,:,:,Rho_]
Can2_RYF = @. RS_YV[:,:,:,By_]^2 / RS_YV[:,:,:,Rho_]
Can2_RZF = @. RS_ZV[:,:,:,Bz_]^2 / RS_ZV[:,:,:,Rho_]

u_LXF = @. LS_XV[:,:,:,Ux_]/LS_XV[:,:,:,Rho_]
u_LYF = @. LS_YV[:,:,:,Uy_]/LS_YV[:,:,:,Rho_]
u_LZF = @. LS_ZV[:,:,:,Uz_]/LS_ZV[:,:,:,Rho_]
u_RXF = @. RS_XV[:,:,:,Ux_]/RS_XV[:,:,:,Rho_]
u_RYF = @. RS_YV[:,:,:,Uy_]/RS_YV[:,:,:,Rho_]
u_RZF = @. RS_ZV[:,:,:,Uz_]/RS_ZV[:,:,:,Rho_]

c_LXF = @. sqrt( 0.5*(Cs2_LXF + Ca2_LXF +
sqrt((Cs2_LXF + Ca2_LXF)^2 - 4.0*Cs2_LXF*Can2_LXF)) )
c_LYF = @. sqrt( 0.5*(Cs2_LYF + Ca2_LYF +
sqrt((Cs2_LYF + Ca2_LYF)^2 - 4.0*Cs2_LYF*Can2_LYF)) )
c_LZF = @. sqrt( 0.5*(Cs2_LZF + Ca2_LZF +
sqrt((Cs2_LZF + Ca2_LZF)^2 - 4.0*Cs2_LZF*Can2_LZF)) )
c_RXF = @. sqrt( 0.5*(Cs2_RXF + Ca2_RXF +
sqrt((Cs2_RXF + Ca2_RXF)^2 - 4.0*Cs2_RXF*Can2_RXF)) )
c_RYF = @. sqrt( 0.5*(Cs2_RYF + Ca2_RYF +
sqrt((Cs2_RYF + Ca2_RYF)^2 - 4.0*Cs2_RYF*Can2_RYF)) )
c_RZF = @. sqrt( 0.5*(Cs2_RZF + Ca2_RZF +
sqrt((Cs2_RZF + Ca2_RZF)^2 - 4.0*Cs2_RZF*Can2_RZF)) )

sLmax_XF = max(0,u_LXF+c_LXF)
sLmin_XF = min(0,u_LXF-c_LXF)
sRmax_XF = max(0,u_RXF+c_RXF)
sRmin_XF = min(0,u_RXF-c_RXF)
sLmax_YF = max(0,u_LYF+c_LYF)
sLmin_YF = min(0,u_LYF-c_LYF)
sRmax_YF = max(0,u_RYF+c_RYF)
sRmin_YF = min(0,u_RYF-c_RYF)
sLmax_ZF = max(0,u_LZF+c_LZF)
sLmin_ZF = min(0,u_LZF-c_LZF)
sRmax_ZF = max(0,u_RZF+c_RZF)
sRmin_ZF = min(0,u_RZF-c_RZF)

smax_XF = max(sLmax_XF,sRmax_XF)
smin_XF = min(sLmin_XF,sRmin_XF)
smax_YF = max(sLmax_YF,sRmax_YF)
smin_YF = min(sLmin_YF,sRmin_YF)
smax_ZF = max(sLmax_ZF,sRmax_ZF)
smin_ZF = min(sLmin_ZF,sRmin_ZF)


# remove dropdims?
bL = dropdims(sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4)
bR = dropdims(sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4)


# HLLE
@. Flux_XV -=
   0.5*(Cmax_XF + Cmin_XF)/(Cmax_XF - Cmin_XF)*(RFlux_XV - LFlux_XV) -
   Cmax_XF*Cmin_XF/(Cmax_XF - Cmin_XF)*(RState_XV - LState_XV)
@. Flux_YV -=
   0.5*(Cmax_YF + Cmin_YF)/(Cmax_YF - Cmin_YF)*(RFlux_YV - LFlux_YV) -
   Cmax_YF*Cmin_YF/(Cmax_YF - Cmin_YF)*(RState_YV - LState_YV)
@. Flux_ZV -=
   0.5*(Cmax_ZF + Cmin_ZF)/(Cmax_ZF - Cmin_ZF)*(RFlux_ZV - LFlux_ZV) -
   Cmax_ZF*Cmin_ZF/(Cmax_ZF - Cmin_ZF)*(RState_ZV - LState_ZV)


 @. Flux_XV[:,:,:,Rho_:Bz_] -= 0.5*(Cmax_XF + Cmin_XF)/(Cmax_XF - Cmin_XF)*
 (RFlux_XV[:,:,:,Rho_:Bz_] - LFlux_XV[:,:,:,Rho_:Bz_]) -
 Cmax_XF*Cmin_XF/(Cmax_XF - Cmin_XF)*(RState_XV[:,:,:,Rho_:Bz_] - LState_XV[:,:,:,Rho_:Bz_])

 @. Flux_XV[:,:,:,E_] -= 0.5*(Cmax_XF + Cmin_XF)/(Cmax_XF - Cmin_XF)*
 (RFlux_XV[:,:,:,E_] - LFlux_XV[:,:,:,E_]) -
 Cmax_XF*Cmin_XF/(Cmax_XF - Cmin_XF)* (
 (RState_XV[:,:,:,P_] / (γ-1) +
 0.5/RState_XV[:,:,:,Rho_]*$dropdims($sum(RState_XV[:,:,:,U_].^2,dims=4);dims=4) +
 0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4)) -
 (LState_XV[:,:,:,P_] / (γ-1) +
 0.5/LState_XV[:,:,:,Rho_]*$dropdims($sum(LState_XV[:,:,:,U_].^2,dims=4);dims=4) +
 0.5*$dropdims($sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4)))

 @. Flux_YV[:,:,:,Rho_:Bz_] -= 0.5*(Cmax_YF + Cmin_YF)/(Cmax_YF - Cmin_YF)*
 (RFlux_YV[:,:,:,Rho_:Bz_] - LFlux_YV[:,:,:,Rho_:Bz_]) -
 Cmax_YF*Cmin_YF/(Cmax_YF - Cmin_YF)*
 (RState_YV[:,:,:,Rho_:Bz_] - LState_YV[:,:,:,Rho_:Bz_])
 @. Flux_YV[:,:,:,E_] -=
 0.5*(Cmax_YF + Cmin_YF)/(Cmax_YF - Cmin_YF)*
 (RFlux_YV[:,:,:,E_] - LFlux_YV[:,:,:,E_]) -
 Cmax_YF*Cmin_YF/(Cmax_YF - Cmin_YF)* (
 (RState_YV[:,:,:,P_] / (γ-1) +
 0.5/RState_YV[:,:,:,Rho_]*$dropdims($sum(RState_YV[:,:,:,U_].^2,dims=4);dims=4) +
 0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4)) -
 (LState_YV[:,:,:,P_] / (γ-1) +
 0.5/LState_YV[:,:,:,Rho_]*$dropdims($sum(LState_YV[:,:,:,U_].^2,dims=4);dims=4) +
 0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4)))

 @. Flux_ZV[:,:,:,Rho_:Bz_] -= 0.5*(Cmax_ZF + Cmin_ZF)/(Cmax_ZF - Cmin_ZF)*
 (RFlux_ZV[:,:,:,Rho_:Bz_] - LFlux_ZV[:,:,:,Rho_:Bz_]) -
 Cmax_ZF*Cmin_ZF/(Cmax_ZF - Cmin_ZF)*
 (RState_ZV[:,:,:,Rho_:Bz_] - LState_ZV[:,:,:,Rho_:Bz_])
 @. Flux_ZV[:,:,:,E_] -=
 0.5*(Cmax_ZF + Cmin_ZF)/(Cmax_ZF - Cmin_ZF)*
 (RFlux_ZV[:,:,:,E_] - LFlux_ZV[:,:,:,E_]) -
 Cmax_ZF*Cmin_ZF/(Cmax_ZF - Cmin_ZF)* (
 (RState_ZV[:,:,:,P_] / (γ-1) +
 0.5/RState_ZV[:,:,:,Rho_]*$dropdims($sum(RState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
 0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4)) -
 (LState_ZV[:,:,:,P_] / (γ-1) +
 0.5/LState_ZV[:,:,:,Rho_]*$dropdims($sum(LState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
 0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4)))

# source
for k = 1:nK, j = 1:nJ, i = 1:nI
   source_GV[i,j,k,Rho_] = 0.0
end


# This is slower than dot fusion, because it is wrong in indexing!
for k = 1:nK, j = 1:nJ, i = 1:nI
   source_GV[i,j,k,Ux_] = -state_GV[i,j,k,Bx_]*DivB[i+1,j+1,k+1]
   source_GV[i,j,k,Uy_] = -state_GV[i,j,k,By_]*DivB[i+1,j+1,k+1]
   source_GV[i,j,k,Uz_] = -state_GV[i,j,k,Bz_]*DivB[i+1,j+1,k+1]
   source_GV[i,j,k,Bx_] =  state_GV[i,j,k,Ux_]*DivB[i+1,j+1,k+1]
   source_GV[i,j,k,By_] =  state_GV[i,j,k,Uy_]*DivB[i+1,j+1,k+1]
   source_GV[i,j,k,Bz_] =  state_GV[i,j,k,Uz_]*DivB[i+1,j+1,k+1]
end

@. source_GV[:,:,:,U_] = -state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_] *
   DivB[iMin:iMax,jMin:jMax,kMin:kMax]

@. source_GV[:,:,:,B_] =  state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_] *
   DivB[iMin:iMax,jMin:jMax,kMin:kMax]

@. source_GV[:,:,:,P_] =-(γ-1.0)*state_GV[iMin:iMax,jMin:jMax,kMin:kMax,P_]*
   DivU[iMin:iMax,jMin:jMax,kMin:kMax]

@. source_GV[:,:,:,E_] =
   -$dropdims($sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_] *
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4);dims=4) *
   DivB[iMin:iMax,jMin:jMax,kMin:kMax]

ub = dropdims(sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].*
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4);dims=4)
for k = 1:nK, j = 1:nJ, i = 1:nI
   source_GV[i,j,k,E_] = -ub[i,j,k]*DivB[i,j,k]
end

# test for simd usage, this may be faster.
@inbounds for k = kMin:kMax
   for j = jMin:jMax
      @simd for i = iMin:iMax
         source_GV[i-nG,j-nG,k-nG,Ux_] = -state_GV[i,j,k,Bx_]*DivB[i,j,k]
         source_GV[i-nG,j-nG,k-nG,Uy_] = -state_GV[i,j,k,By_]*DivB[i,j,k]
         source_GV[i-nG,j-nG,k-nG,Uz_] = -state_GV[i,j,k,Bz_]*DivB[i,j,k]
         source_GV[i-nG,j-nG,k-nG,Bx_] =  state_GV[i,j,k,Ux_]*DivB[i,j,k]
         source_GV[i-nG,j-nG,k-nG,By_] =  state_GV[i,j,k,Uy_]*DivB[i,j,k]
         source_GV[i-nG,j-nG,k-nG,Bz_] =  state_GV[i,j,k,Uz_]*DivB[i,j,k]
      end
   end
end

function calc_source!(param::Param, state_GV::Array{Float64,4},
   source_GV::Array{Float64,4})

   nVar, nG = param.nVar, param.nG
   x, y, z = param.x, param.y, param.z

   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax

   # I can preallocate DivB and DivU!
   nI,nJ,nK = param.nI, param.nJ, param.nK


   # Calculate divergence of B using central difference
   DivB = divergence_ndgrid(x,y,z,state_GV[:,:,:,B_])

   source_GV[:,:,:,Rho_] .= 0.0

   @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
      source_GV[i-nG,j-nG,k-nG,Ux_] = -state_GV[i,j,k,Bx_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Uy_] = -state_GV[i,j,k,By_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Uz_] = -state_GV[i,j,k,Bz_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Bx_] =  state_GV[i,j,k,Ux_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,By_] =  state_GV[i,j,k,Uy_]*DivB[i,j,k]
      source_GV[i-nG,j-nG,k-nG,Bz_] =  state_GV[i,j,k,Uz_]*DivB[i,j,k]
   end

   # Calculate divergence of U
   DivU = divergence_ndgrid(x,y,z,state_GV[:,:,:,U_]./state_GV[:,:,:,Rho_])

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

# Update state
@. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,:] =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax,:] - dt*(
   (Flux_XV[2:end,:,:,:] - Flux_XV[1:end-1,:,:,:])/CellSize_D[1] +
   (Flux_YV[:,2:end,:,:] - Flux_YV[:,1:end-1,:,:])/CellSize_D[2] +
   (Flux_ZV[:,:,2:end,:] - Flux_ZV[:,:,1:end-1,:])/CellSize_D[3] -
   source_GV)


@. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_:Bz_] =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_:Bz_] - dt*(
   (Flux_XV[2:end,:,:,Rho_:Bz_] - Flux_XV[1:end-1,:,:,Rho_:Bz_])
   /CellSize_D[1] +
   (Flux_YV[:,2:end,:,Rho_:Bz_] - Flux_YV[:,1:end-1,:,Rho_:Bz_])
   /CellSize_D[2] +
   (Flux_ZV[:,:,2:end,Rho_:Bz_] - Flux_ZV[:,:,1:end-1,Rho_:Bz_])
   /CellSize_D[3] +
   source_GV[:,:,:,Rho_:Bz_])

@. state_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax,P_]/(γ-1) +
   0.5/state_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_]*
   $dropdims($sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,dims=4)
   ;dims=4) + 0.5*
   $dropdims($sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4);
   dims=4)

@. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] =
   state_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] - dt*(
   (Flux_XV[2:end,:,:,E_] - Flux_XV[1:end-1,:,:,E_])/CellSize_D[1] +
   (Flux_YV[:,2:end,:,E_] - Flux_YV[:,1:end-1,:,E_])/CellSize_D[2] +
   (Flux_ZV[:,:,2:end,E_] - Flux_ZV[:,:,1:end-1,E_])/CellSize_D[3] +
   source_GV[:,:,:,E_])

@. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,P_] = (γ-1.0)*
   (stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] -
   0.5/stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_]*
   $dropdims($sum(stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,
   dims=4);dims=4) - 0.5*
   $dropdims($sum(stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4)
   ;dims=4) )


u = sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,dims=4)
b = sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4) # takes time

u = sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,dims=4)
b = sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4)

   function update_state!(param::Param, state_GV::Array{Float64,4}, dt::Float64,
      faceFlux::FaceFlux, source_GV::Array{Float64,4})

      # In local timestepping, dt can be an array!

      Flux_XV = faceFlux.Flux_XV
      Flux_YV = faceFlux.Flux_YV
      Flux_ZV = faceFlux.Flux_ZV

      #stateNew_GV = similar(state_GV)
      stateNew_GV = copy(state_GV)

      CellSize_D = param.CellSize_D

      iMin, iMax, jMin, jMax, kMin, kMax =
      param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax
      nVar = param.nVar
      nI,nJ,nK = param.nI, param.nJ, param.nK

      if param.TypeGrid == "Cartesian"
         # No need for volume and face if the grid is uniform Cartesian
         if !param.UseConservative
            for iVar=1:nVar, k=1:nK, j=1:nJ, i=1:nI
               stateNew_GV[i+1,j+1,k+1,iVar] = state_GV[i+1,j+1,k+1,iVar] - dt*(
               (Flux_XV[i+1,j,k,iVar] - Flux_XV[i,j,k,iVar])/CellSize_D[1] +
               (Flux_YV[i,j+1,k,iVar] - Flux_YV[i,j,k,iVar])/CellSize_D[2] +
               (Flux_ZV[i,j,k+1,iVar] - Flux_ZV[i,j,k,iVar])/CellSize_D[3] +
               source_GV[i,j,k,iVar])
            end
         else
            for iVar=Rho_:Bz_, k=1:nK, j=1:nJ, i=1:nI
               stateNew_GV[i+1,j+1,k+1,iVar] = state_GV[i+1,j+1,k+1,iVar] - dt*(
               (Flux_XV[i+1,j,k,iVar] - Flux_XV[i,j,k,iVar])/CellSize_D[1] +
               (Flux_YV[i,j+1,k,iVar] - Flux_YV[i,j,k,iVar])/CellSize_D[2] +
               (Flux_ZV[i,j,k+1,iVar] - Flux_ZV[i,j,k+1,iVar])/CellSize_D[3] +
               source_GV[i,j,k,iVar])
            end

            u = dropdims(sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,dims=4);dims=4)
            b = dropdims(sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4);dims=4)
            for k=1:nK, j=1:nJ, i=1:nI
               state_GV[i+1,j+1,k+1,E_] = state_GV[i+1,j+1,k+1,P_]/(γ-1.0) +
                  0.5/state_GV[i+1,j+1,k+1,Rho_]*u[i,j,k] + 0.5*b[i,j,k]
            end

            for k=1:nK, j=1:nJ, i=1:nI
               stateNew_GV[i+1,j+1,k+1,E_] = state_GV[i+1,j+1,k+1,E_] - dt*(
               (Flux_XV[i+1,j,k,E_] - Flux_XV[i,j,k,E_])/CellSize_D[1] +
               (Flux_YV[i,j+1,k,E_] - Flux_YV[i,j,k,E_])/CellSize_D[2] +
               (Flux_ZV[i,j,k+1,E_] - Flux_ZV[i,j,k,E_])/CellSize_D[3] +
               source_GV[i,j,k,E_])
            end

            u = dropdims(sum(stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,dims=4);dims=4)
            b = dropdims(sum(stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4);dims=4)

            for k=1:nK, j=1:nJ, i=1:nI
               stateNew_GV[i+1,j+1,k+1,P_] = (γ-1.0)*(stateNew_GV[i+1,j+1,k+1,E_] -
                  0.5/stateNew_GV[i+1,j+1,k+1,Rho_]*u[i,j,k] - 0.5*b[i,j,k])
            end
         end
      else
         # Need volume and face
         stateNew_GV .= 0.0
      end

      state_GV .= stateNew_GV

   end

   function update_state!(param::Param, state_GV::Array{Float64,4}, dt::Float64,
      faceFlux::FaceFlux, source_GV::Array{Float64,4})

      # In local timestepping, dt can be an array!

      Flux_XV = faceFlux.Flux_XV
      Flux_YV = faceFlux.Flux_YV
      Flux_ZV = faceFlux.Flux_ZV

      #stateNew_GV = similar(state_GV)
      stateNew_GV = copy(state_GV)

      CellSize_D = param.CellSize_D

      iMin, iMax, jMin, jMax, kMin, kMax =
      param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax
      nVar = param.nVar
      nI,nJ,nK = param.nI, param.nJ, param.nK

      if param.TypeGrid == "Cartesian"
         # No need for volume and face if the grid is uniform Cartesian
         if !param.UseConservative
            @. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,:] =
               state_GV[iMin:iMax,jMin:jMax,kMin:kMax,:] - dt*(
               (Flux_XV[2:end,:,:,:] - Flux_XV[1:end-1,:,:,:])/CellSize_D[1] +
               (Flux_YV[:,2:end,:,:] - Flux_YV[:,1:end-1,:,:])/CellSize_D[2] +
               (Flux_ZV[:,:,2:end,:] - Flux_ZV[:,:,1:end-1,:])/CellSize_D[3] +
               source_GV)
         else
            @. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_:Bz_] =
               state_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_:Bz_] - dt*(
               (Flux_XV[2:end,:,:,Rho_:Bz_] - Flux_XV[1:end-1,:,:,Rho_:Bz_])
               /CellSize_D[1] +
               (Flux_YV[:,2:end,:,Rho_:Bz_] - Flux_YV[:,1:end-1,:,Rho_:Bz_])
               /CellSize_D[2] +
               (Flux_ZV[:,:,2:end,Rho_:Bz_] - Flux_ZV[:,:,1:end-1,Rho_:Bz_])
               /CellSize_D[3] +
               source_GV[:,:,:,Rho_:Bz_])

            @. state_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] =
               state_GV[iMin:iMax,jMin:jMax,kMin:kMax,P_]/(γ-1) +
               0.5/state_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_]*
               $dropdims($sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,dims=4)
               ;dims=4) + 0.5*
               $dropdims($sum(state_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4);
               dims=4)

            @. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] =
               state_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] - dt*(
               (Flux_XV[2:end,:,:,E_] - Flux_XV[1:end-1,:,:,E_])/CellSize_D[1] +
               (Flux_YV[:,2:end,:,E_] - Flux_YV[:,1:end-1,:,E_])/CellSize_D[2] +
               (Flux_ZV[:,:,2:end,E_] - Flux_ZV[:,:,1:end-1,E_])/CellSize_D[3] +
               source_GV[:,:,:,E_])

            @. stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,P_] = (γ-1.0)*
               (stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,E_] -
               0.5/stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,Rho_]*
               $dropdims($sum(stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,U_].^2,
               dims=4);dims=4) - 0.5*
               $dropdims($sum(stateNew_GV[iMin:iMax,jMin:jMax,kMin:kMax,B_],dims=4)
               ;dims=4) )
         end
      else
         # Need volume and face
         stateNew_GV .= 0.0
      end

      state_GV .= stateNew_GV

   end

   time_G .= CFL ./ (
   (nI > 1).*Cmax_XF./CellSize_D[1] .+
   (nJ > 1).*Cmax_YF./CellSize_D[2] .+
   (nK > 1).*Cmax_ZF./CellSize_D[3])

# Divergence
@. px[2:n-1,:,:] = (var[3:n,:,:,1] - var[1:n-2,:,:,1])/(hx[3:n] - hx[1:n-2])
@. qy[:,2:n-1,:] = (var[:,3:n,:,2] - var[:,1:n-2,:,2])/(hy[3:n] - hy[1:n-2])
@. rz[:,:,2:n-1] = (var[:,:,3:n,3] - var[:,:,1:n-2,3])/(hz[3:n] - hz[1:n-2])

# Boundary

# periodic
state_GV[1:nG,jMin:jMax,kMin:kMax,:] .= state_GV[iMax-nG+1:iMax,jMin:jMax,kMin:kMax,:]
state_GV[iMax+1:iMaxAll,jMin:jMax,kMin:kMax,:] .= state_GV[iMin:iMin+nG-1,jMin:jMax,kMin:kMax,:]
state_GV[iMin:iMax,1:nG,kMin:kMax,:] .= state_GV[iMin:iMax,jMax-nG+1:jMax,kMin:kMax,:]
state_GV[iMin:iMax,jMax+1:jMaxAll,kMin:kMax,:] .= state_GV[iMin:iMax,jMin:jMin+nG-1,kMin:kMax,:]
state_GV[iMin:iMax,jMin:jMax,1:nG,:] .= state_GV[iMin:iMax,jMin:jMax,kMax-nG+1:kMax,:]
state_GV[iMin:iMax,jMin:jMax,kMax+1:kMaxAll,:] .= state_GV[iMin:iMax,jMin:jMax,kMin:kMin+nG-1,:]
# float
state_GV[1:nG,jMin:jMax,kMin:kMax,:] .=
   state_GV[iMin:iMin,jMin:jMax,kMin:kMax,:]
state_GV[iMax+1:iMaxAll,jMin:jMax,kMin:kMax,:] .=
   state_GV[iMax:iMax,jMin:jMax,kMin:kMax,:]

state_GV[iMin:iMax,1:nG,kMin:kMax,:] .=
   state_GV[iMin:iMax,jMin:jMin,kMin:kMax,:]
state_GV[iMin:iMax,jMax+1:jMaxAll,kMin:kMax,:] .=
   state_GV[iMin:iMax,jMax:jMax,kMin:kMax,:]

state_GV[iMin:iMax,jMin:jMax,1:nG,:] .=
   state_GV[iMin:iMax,jMin:jMax,kMin:kMin,:]
state_GV[iMin:iMax,jMin:jMax,kMax+1:kMaxAll,:] .=
   state_GV[iMin:iMax,jMin:jMax,kMax:kMax,:]

# calc_timestep
Cmax_XG = max.(Cmax_XF[2:end,:,:],Cmax_XF[1:end-1,:,:]) # these allocations are not needed!
Cmax_YG = max.(Cmax_YF[:,2:end,:],Cmax_YF[:,1:end-1,:])
Cmax_ZG = max.(Cmax_ZF[:,:,2:end],Cmax_ZF[:,:,1:end-1])

time_G = CFL ./ (
   (nI > 1).*Cmax_XG./CellSize_D[1] .+
   (nJ > 1).*Cmax_YG./CellSize_D[2] .+
   (nK > 1).*Cmax_ZG./CellSize_D[3])

   # flux module using dot fusion.
   module Flux

   export calc_face_flux!, init_flux, SpeedFlux, FaceFlux, FaceFluxLR

   using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_
   using ..Parameters: γ
   using ..FaceValue: FaceState

   struct SpeedFlux
      Cmax_XF::Array{Float64,3}
      Cmax_YF::Array{Float64,3}
      Cmax_ZF::Array{Float64,3}
   end

   struct FaceFlux
      Flux_XV::Array{Float64,4}
      Flux_YV::Array{Float64,4}
      Flux_ZV::Array{Float64,4}
   end

   struct FaceFluxLR
      LFlux_XV::Array{Float64,4}
      RFlux_XV::Array{Float64,4}
      LFlux_YV::Array{Float64,4}
      RFlux_YV::Array{Float64,4}
      LFlux_ZV::Array{Float64,4}
      RFlux_ZV::Array{Float64,4}
   end

   function init_flux(param::Param)

      GridSize = param.GridSize
      nVar = param.nVar

      Flux_XV = Array{Float64,4}(undef,GridSize+[1,0,0]...,nVar)
      Flux_YV = Array{Float64,4}(undef,GridSize+[0,1,0]...,nVar)
      Flux_ZV = Array{Float64,4}(undef,GridSize+[0,0,1]...,nVar)

      LFlux_XV = Array{Float64,4}(undef,GridSize+[1,0,0]...,nVar)
      RFlux_XV = Array{Float64,4}(undef,GridSize+[1,0,0]...,nVar)
      LFlux_YV = Array{Float64,4}(undef,GridSize+[0,1,0]...,nVar)
      RFlux_YV = Array{Float64,4}(undef,GridSize+[0,1,0]...,nVar)
      LFlux_ZV = Array{Float64,4}(undef,GridSize+[0,0,1]...,nVar)
      RFlux_ZV = Array{Float64,4}(undef,GridSize+[0,0,1]...,nVar)

      Cmax_XF = Array{Float64,3}(undef,GridSize+[1,0,0]...)
      Cmax_YF = Array{Float64,3}(undef,GridSize+[0,1,0]...)
      Cmax_ZF = Array{Float64,3}(undef,GridSize+[0,0,1]...)

      faceFluxLR = FaceFluxLR(LFlux_XV,RFlux_XV,LFlux_YV,RFlux_YV,LFlux_ZV, RFlux_ZV)
      faceFlux = FaceFlux(Flux_XV, Flux_YV, Flux_ZV)
      speedFlux = SpeedFlux(Cmax_XF, Cmax_YF, Cmax_ZF)

      return faceFluxLR, faceFlux, speedFlux
   end

   function calc_face_flux!(param::Param, faceValue::FaceState,
      faceFlux::FaceFlux, speedFlux::SpeedFlux, faceFluxLR::FaceFluxLR)

      nVar = param.nVar
      GridSize = param.GridSize

      if param.Scheme == "Rusanov"
         get_physical_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)
      elseif param.Scheme == "HLLE"
         get_HLLE_flux(param, faceValue, faceFlux, speedFlux)
      end

      return
   end

   function get_physical_flux!(param::Param, faceValue::FaceState,
      faceFlux::FaceFlux, speedFlux::SpeedFlux, faceFluxLR::FaceFluxLR)

      LState_XV = faceValue.LState_XV
      RState_XV = faceValue.RState_XV
      LState_YV = faceValue.LState_YV
      RState_YV = faceValue.RState_YV
      LState_ZV = faceValue.LState_ZV
      RState_ZV = faceValue.RState_ZV

      LFlux_XV = faceFluxLR.LFlux_XV
      RFlux_XV = faceFluxLR.RFlux_XV
      LFlux_YV = faceFluxLR.LFlux_YV
      RFlux_YV = faceFluxLR.RFlux_YV
      LFlux_ZV = faceFluxLR.LFlux_ZV
      RFlux_ZV = faceFluxLR.RFlux_ZV

      Flux_XV = faceFlux.Flux_XV
      Flux_YV = faceFlux.Flux_YV
      Flux_ZV = faceFlux.Flux_ZV

      GridSize = param.GridSize
      nVar = param.nVar

      nI,nJ,nK = param.nI, param.nJ, param.nK

      # Density flux
      @. LFlux_XV[:,:,:,Rho_] = LState_XV[:,:,:,Ux_]
      @. RFlux_XV[:,:,:,Rho_] = RState_XV[:,:,:,Ux_]
      @. LFlux_YV[:,:,:,Rho_] = LState_YV[:,:,:,Uy_]
      @. RFlux_YV[:,:,:,Rho_] = RState_YV[:,:,:,Uy_]
      @. LFlux_ZV[:,:,:,Rho_] = LState_ZV[:,:,:,Uz_]
      @. RFlux_ZV[:,:,:,Rho_] = RState_ZV[:,:,:,Uz_]

      # Momentum flux
      @. LFlux_XV[:,:,:,Ux_] = LState_XV[:,:,:,Ux_]^2 / LState_XV[:,:,:,Rho_] + LState_XV[:,:,:,P_] + 0.5*$dropdims($sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4) - LState_XV[:,:,:,Bx_]^2
      @. RFlux_XV[:,:,:,Ux_] = RState_XV[:,:,:,Ux_]^2 / RState_XV[:,:,:,Rho_] + RState_XV[:,:,:,P_] + 0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4) - RState_XV[:,:,:,Bx_]^2
      @. LFlux_XV[:,:,:,Uy_] = LState_XV[:,:,:,Ux_] * LState_XV[:,:,:,Uy_] / LState_XV[:,:,:,Rho_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,By_]
      @. RFlux_XV[:,:,:,Uy_] = RState_XV[:,:,:,Ux_] * RState_XV[:,:,:,Uy_] / RState_XV[:,:,:,Rho_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,By_]
      @. LFlux_XV[:,:,:,Uz_] = LState_XV[:,:,:,Ux_] * LState_XV[:,:,:,Uz_] / LState_XV[:,:,:,Rho_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,Bz_]
      @. RFlux_XV[:,:,:,Uz_] = RState_XV[:,:,:,Ux_] * RState_XV[:,:,:,Uz_] / RState_XV[:,:,:,Rho_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,Bz_]

      @. LFlux_YV[:,:,:,Ux_] = LState_YV[:,:,:,Uy_] * LState_YV[:,:,:,Ux_] / LState_YV[:,:,:,Rho_] - LState_YV[:,:,:,By_]*LState_YV[:,:,:,Bx_]
      @. RFlux_YV[:,:,:,Ux_] = RState_YV[:,:,:,Uy_] * RState_YV[:,:,:,Ux_] / RState_YV[:,:,:,Rho_] - RState_YV[:,:,:,By_]*RState_YV[:,:,:,Bx_]
      @. LFlux_YV[:,:,:,Uy_] = LState_YV[:,:,:,Uy_]^2 / LState_YV[:,:,:,Rho_] + LState_YV[:,:,:,P_] + 0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4) - LState_YV[:,:,:,By_]^2
      @. RFlux_YV[:,:,:,Uy_] = RState_YV[:,:,:,Uy_]^2 / RState_YV[:,:,:,Rho_] + RState_YV[:,:,:,P_] + 0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4) - RState_YV[:,:,:,By_]^2
      @. LFlux_YV[:,:,:,Uz_] = LState_YV[:,:,:,Uy_] * LState_YV[:,:,:,Uz_] / LState_YV[:,:,:,Rho_] - LState_YV[:,:,:,Bx_]*LState_YV[:,:,:,Bz_]
      @. RFlux_YV[:,:,:,Uz_] = RState_YV[:,:,:,Uy_] * RState_YV[:,:,:,Uz_] / RState_YV[:,:,:,Rho_] - RState_YV[:,:,:,Bx_]*RState_YV[:,:,:,Bz_]

      @. LFlux_ZV[:,:,:,Ux_] = LState_ZV[:,:,:,Uz_] * LState_ZV[:,:,:,Ux_] / LState_ZV[:,:,:,Rho_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,Bx_]
      @. RFlux_ZV[:,:,:,Ux_] = RState_ZV[:,:,:,Uz_] * RState_ZV[:,:,:,Ux_] / RState_ZV[:,:,:,Rho_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,Bx_]
      @. LFlux_ZV[:,:,:,Uy_] = LState_ZV[:,:,:,Uz_] * LState_ZV[:,:,:,Uy_] / LState_ZV[:,:,:,Rho_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,By_]
      @. RFlux_ZV[:,:,:,Uy_] = RState_ZV[:,:,:,Uz_] * RState_ZV[:,:,:,Uy_] / RState_ZV[:,:,:,Rho_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,By_]
      @. LFlux_ZV[:,:,:,Uz_] = LState_ZV[:,:,:,Uz_]^2 / LState_ZV[:,:,:,Rho_] + LState_ZV[:,:,:,P_] + 0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4) - LState_ZV[:,:,:,Bz_]^2
      @. RFlux_ZV[:,:,:,Uz_] = RState_ZV[:,:,:,Uz_]^2 / RState_ZV[:,:,:,Rho_] + RState_ZV[:,:,:,P_] + 0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4) - RState_ZV[:,:,:,Bz_]^2

      # Magnetic flux
      @. LFlux_XV[:,:,:,Bx_] = 0.0
      @. RFlux_XV[:,:,:,Bx_] = 0.0
      @. LFlux_XV[:,:,:,By_] = LState_XV[:,:,:,Ux_]*LState_XV[:,:,:,By_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,Uy_]
      @. RFlux_XV[:,:,:,By_] = RState_XV[:,:,:,Ux_]*RState_XV[:,:,:,By_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,Uy_]
      @. LFlux_XV[:,:,:,Bz_] = LState_XV[:,:,:,Ux_]*LState_XV[:,:,:,Bz_] - LState_XV[:,:,:,Bx_]*LState_XV[:,:,:,Uz_]
      @. RFlux_XV[:,:,:,Bz_] = RState_XV[:,:,:,Ux_]*RState_XV[:,:,:,Bz_] - RState_XV[:,:,:,Bx_]*RState_XV[:,:,:,Uz_]

      @. LFlux_YV[:,:,:,Bx_] = LState_YV[:,:,:,Uy_]*LState_YV[:,:,:,Bx_] - LState_YV[:,:,:,By_]*LState_YV[:,:,:,Ux_]
      @. RFlux_YV[:,:,:,Bx_] = RState_YV[:,:,:,Uy_]*RState_YV[:,:,:,Bx_] - RState_YV[:,:,:,By_]*RState_YV[:,:,:,Ux_]
      @. LFlux_YV[:,:,:,By_] = 0.0
      @. RFlux_YV[:,:,:,By_] = 0.0
      @. LFlux_YV[:,:,:,Bz_] = LState_YV[:,:,:,Uy_]*LState_YV[:,:,:,Bz_] - LState_YV[:,:,:,By_]*LState_YV[:,:,:,Uz_]
      @. RFlux_YV[:,:,:,Bz_] = RState_YV[:,:,:,Uy_]*RState_YV[:,:,:,Bz_] - RState_YV[:,:,:,By_]*RState_YV[:,:,:,Uz_]

      @. LFlux_ZV[:,:,:,Bx_] = LState_ZV[:,:,:,Uz_]*LState_ZV[:,:,:,Bx_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,Ux_]
      @. RFlux_ZV[:,:,:,Bx_] = RState_ZV[:,:,:,Uz_]*RState_ZV[:,:,:,Bx_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,Ux_]
      @. LFlux_ZV[:,:,:,By_] = LState_ZV[:,:,:,Uz_]*LState_ZV[:,:,:,By_] - LState_ZV[:,:,:,Bz_]*LState_ZV[:,:,:,Uy_]
      @. RFlux_ZV[:,:,:,By_] = RState_ZV[:,:,:,Uz_]*RState_ZV[:,:,:,By_] - RState_ZV[:,:,:,Bz_]*RState_ZV[:,:,:,Uy_]
      @. LFlux_ZV[:,:,:,Bz_] = 0.0
      @. RFlux_ZV[:,:,:,Bz_] = 0.0

      # Pressure flux / energy flux
      if !param.UseConservative
         for k = 1:nK+1, j = 1:nJ, i = 1:nI
            @. LFlux_XV[:,:,:,P_] = LState_XV[:,:,:,P_]*LState_XV[:,:,:,Ux_]/LState_XV[:,:,:,Rho_]
            @. RFlux_XV[:,:,:,P_] = RState_XV[:,:,:,P_]*RState_XV[:,:,:,Ux_]/RState_XV[:,:,:,Rho_]
            @. LFlux_YV[:,:,:,P_] = LState_YV[:,:,:,P_]*LState_YV[:,:,:,Uy_]/LState_YV[:,:,:,Rho_]
            @. RFlux_YV[:,:,:,P_] = RState_YV[:,:,:,P_]*RState_YV[:,:,:,Uy_]/RState_YV[:,:,:,Rho_]
            @. LFlux_ZV[:,:,:,P_] = LState_ZV[:,:,:,P_]*LState_ZV[:,:,:,Uz_]/LState_ZV[:,:,:,Rho_]
            @. RFlux_ZV[:,:,:,P_] = RState_ZV[:,:,:,P_]*RState_ZV[:,:,:,Uz_]/RState_ZV[:,:,:,Rho_]
         end
      else
         # Currently use the same index for pressure/energy
         @. LFlux_XV[:,:,:,E_] = LState_XV[:,:,:,Ux_]/LState_XV[:,:,:,Rho_]*
               ((LState_XV[:,:,:,P_]/(γ-1) + 0.5/LState_XV[:,:,:,Rho_]*
               $dropdims($sum(LState_XV[:,:,:,U_].^2,dims=4);dims=4) +
               0.5*$dropdims($sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4) + LState_XV[:,:,:,P_] +
               0.5*$dropdims($sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4))) -
               $dropdims($sum(LState_XV[:,:,:,U_].*LState_XV[:,:,:,B_],dims=4);dims=4)*
               LState_XV[:,:,:,Bx_]

         @. RFlux_XV[:,:,:,E_] = RState_XV[:,:,:,Ux_]/RState_XV[:,:,:,Rho_]*
               ((RState_XV[:,:,:,P_]/(γ-1) + 0.5/RState_XV[:,:,:,Rho_]*
               $dropdims($sum(RState_XV[:,:,:,U_].^2,dims=4);dims=4) +
               0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4) + RState_XV[:,:,:,P_] +
               0.5*$dropdims($sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4))) -
               $dropdims($sum(RState_XV[:,:,:,U_].*RState_XV[:,:,:,B_],dims=4);dims=4)*
               RState_XV[:,:,:,Bx_]

         @. LFlux_YV[:,:,:,E_] = LState_YV[:,:,:,Uy_]/LState_YV[:,:,:,Rho_]*
               ((LState_YV[:,:,:,P_]/(γ-1) + 0.5/LState_YV[:,:,:,Rho_]*
               $dropdims($sum(LState_YV[:,:,:,U_].^2,dims=4);dims=4) +
               0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4) + LState_YV[:,:,:,P_] +
               0.5*$dropdims($sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4))) -
               $dropdims($sum(LState_YV[:,:,:,U_].*LState_YV[:,:,:,B_],dims=4);dims=4)*
               LState_YV[:,:,:,By_]

         @. RFlux_YV[:,:,:,E_] = RState_YV[:,:,:,Uy_]/RState_YV[:,:,:,Rho_]*
               ((RState_YV[:,:,:,P_]/(γ-1) +
               0.5/RState_YV[:,:,:,Rho_]*$dropdims($sum(RState_YV[:,:,:,U_].^2,dims=4);dims=4) +
               0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4) + RState_YV[:,:,:,P_] +
               0.5*$dropdims($sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4))) -
               $dropdims($sum(RState_YV[:,:,:,U_].*RState_YV[:,:,:,B_],dims=4);dims=4)*
               RState_YV[:,:,:,By_]

         @. LFlux_ZV[:,:,:,E_] = LState_ZV[:,:,:,Uz_]/LState_ZV[:,:,:,Rho_]*
               ((LState_ZV[:,:,:,P_]/(γ-1) + 0.5/LState_ZV[:,:,:,Rho_]*
               $dropdims($sum(LState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
               0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4) + LState_ZV[:,:,:,P_] +
               0.5*$dropdims($sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4))) -
               $dropdims($sum(LState_ZV[:,:,:,U_].*LState_ZV[:,:,:,B_],dims=4);dims=4)*
               LState_ZV[:,:,:,Bz_]

         @. RFlux_ZV[:,:,:,E_] = RState_ZV[:,:,:,Uz_]/RState_ZV[:,:,:,Rho_]*
            ((RState_ZV[:,:,:,P_]/(γ-1) + 0.5/RState_ZV[:,:,:,Rho_]*
            $dropdims($sum(RState_ZV[:,:,:,U_].^2,dims=4);dims=4) +
            0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4) + RState_ZV[:,:,:,P_] +
            0.5*$dropdims($sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4))) -
            $dropdims($sum(RState_ZV[:,:,:,U_].*RState_ZV[:,:,:,B_],dims=4);dims=4)*
            RState_ZV[:,:,:,Bz_]
      end

      # Collect all the physical fluxes
      @. Flux_XV = 0.5 * (LFlux_XV + RFlux_XV)
      @. Flux_YV = 0.5 * (LFlux_YV + RFlux_YV)
      @. Flux_ZV = 0.5 * (LFlux_ZV + RFlux_ZV)

      # Cmax_XF should be passed in?
      add_numerical_flux!(param, faceValue, faceFlux, speedFlux)
   end

   function add_numerical_flux!(param::Param, faceValue::FaceState,
      faceFlux::FaceFlux, speedFlux::SpeedFlux)

      nI,nJ,nK = param.nI, param.nJ, param.nK
      nVar = param.nVar

      LState_XV = faceValue.LState_XV
      RState_XV = faceValue.RState_XV
      LState_YV = faceValue.LState_YV
      RState_YV = faceValue.RState_YV
      LState_ZV = faceValue.LState_ZV
      RState_ZV = faceValue.RState_ZV

      Flux_XV = faceFlux.Flux_XV
      Flux_YV = faceFlux.Flux_YV
      Flux_ZV = faceFlux.Flux_ZV

      Cmax_XF = speedFlux.Cmax_XF
      Cmax_YF = speedFlux.Cmax_YF
      Cmax_ZF = speedFlux.Cmax_ZF

      if param.Scheme == "Rusanov"
         get_speed_max!(param, faceValue, speedFlux)

         if !param.UseConservative
            for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
               Flux_XV[i,j,k,iVar] -= 0.5*Cmax_XF[i,j,k]*(RState_XV[i,j,k,iVar] -
                  LState_XV[i,j,k,iVar])
               Flux_YV[i,j,k,iVar] -= 0.5*Cmax_YF[i,j,k]*(RState_YV[i,j,k,iVar] -
                  LState_YV[i,j,k,iVar])
               Flux_ZV[i,j,k,iVar] -= 0.5*Cmax_ZF[i,j,k]*(RState_ZV[i,j,k,iVar] -
                  LState_ZV[i,j,k,iVar])
            end
         else
            # If I solve energy equation instead of pressure, there's
            # duplicate calculation above, even though the expression
            # looks compact. That's why I use an if-else statement
            for iVar = Rho_:Bz_, k = 1:nK, j = 1:nJ, i = 1:nI+1
               Flux_XV[i,j,k,iVar] -= 0.5*Cmax_XF[i,j,k]*(RState_XV[i,j,k,iVar] -
                  LState_XV[i,j,k,iVar])
            end

            uL = dropdims(sum(LState_XV[:,:,:,U_].^2,dims=4);dims=4)
            bL = dropdims(sum(LState_XV[:,:,:,B_].^2,dims=4);dims=4)
            uR = dropdims(sum(RState_XV[:,:,:,U_].^2,dims=4);dims=4)
            bR = dropdims(sum(RState_XV[:,:,:,B_].^2,dims=4);dims=4)

            for k = 1:nK, j = 1:nJ, i = 1:nI+1
               Flux_XV[i,j,k,E_] -= 0.5*Cmax_XF[i,j,k]*(
                  (RState_XV[i,j,k,P_]/(γ-1) + 0.5/RState_XV[i,j,k,Rho_]*uR[i,j,k]
                  + 0.5*bR[i,j,k]) -
                  (LState_XV[i,j,k,P_]/(γ-1) + 0.5/LState_XV[i,j,k,Rho_]*uL[i,j,k]
                  + 0.5*bL[i,j,k]))
            end

            for iVar = Rho_:Bz_, k = 1:nK, j = 1:nJ+1, i = 1:nI
               Flux_YV[i,j,k,iVar] -= 0.5*Cmax_YF[i,j,k]*(RState_YV[i,j,k,iVar] -
                  LState_YV[i,j,k,iVar])
            end

            uL = dropdims(sum(LState_YV[:,:,:,U_].^2,dims=4);dims=4)
            bL = dropdims(sum(LState_YV[:,:,:,B_].^2,dims=4);dims=4)
            uR = dropdims(sum(RState_YV[:,:,:,U_].^2,dims=4);dims=4)
            bR = dropdims(sum(RState_YV[:,:,:,B_].^2,dims=4);dims=4)

            for k = 1:nK, j = 1:nJ+1, i = 1:nI
               Flux_YV[i,j,k,E_] -= 0.5*Cmax_YF[i,j,k]*(
                  (RState_YV[i,j,k,P_]/(γ-1) + 0.5/RState_YV[i,j,k,Rho_]*uR[i,j,k]
                  + 0.5*bR[i,j,k]) -
                  (LState_YV[i,j,k,P_]/(γ-1) + 0.5/LState_YV[i,j,k,Rho_]*uL[i,j,k]
                  + 0.5*bL[i,j,k]))
            end

            for iVar = Rho_:Bz_, k = 1:nK+1, j = 1:nJ, i = 1:nI
               Flux_ZV[i,j,k,iVar] -= 0.5*Cmax_ZF[i,j,k]*(RState_ZV[i,j,k,iVar] -
                  LState_ZV[i,j,k,iVar])
            end

            uL = dropdims(sum(LState_ZV[:,:,:,U_].^2,dims=4);dims=4)
            bL = dropdims(sum(LState_ZV[:,:,:,B_].^2,dims=4);dims=4)
            uR = dropdims(sum(RState_ZV[:,:,:,U_].^2,dims=4);dims=4)
            bR = dropdims(sum(RState_ZV[:,:,:,B_].^2,dims=4);dims=4)

            for k = 1:nK+1, j = 1:nJ, i = 1:nI
               Flux_ZV[i,j,k,E_] -= 0.5*Cmax_ZF[i,j,k]*(
                  (RState_ZV[i,j,k,P_]/(γ-1) + 0.5/RState_ZV[i,j,k,Rho_]*uR[i,j,k]
                  + 0.5*bR[i,j,k]) -
                  (LState_ZV[i,j,k,P_]/(γ-1) + 0.5/LState_ZV[i,j,k,Rho_]*uL[i,j,k]
                  + 0.5*bL[i,j,k]))
            end

         end
      end
   end

   function get_HLLE_flux()

   end

   ```Calculate the maximum speed in each direction.```
   function get_speed_max!(param::Param,faceValue::FaceState,speedFlux::SpeedFlux)
      nI, nJ, nK = param.nI, param.nJ, param.nK
      # Aliases
      LS_XV, RS_XV = faceValue.LState_XV, faceValue.RState_XV
      LS_YV, RS_YV = faceValue.LState_YV, faceValue.RState_YV
      LS_ZV, RS_ZV = faceValue.LState_ZV, faceValue.RState_ZV

      Cmax_XF = speedFlux.Cmax_XF
      Cmax_YF = speedFlux.Cmax_YF
      Cmax_ZF = speedFlux.Cmax_ZF

      for k = 1:nK, j = 1:nJ, i = 1:nI+1
         Cs2_XF = γ*(LS_XV[i,j,k,P_] + RS_XV[i,j,k,P_]) /
            (LS_XV[i,j,k,Rho_] + RS_XV[i,j,k,Rho_])
         Ca2_XF = ( (LS_XV[i,j,k,Bx_] + RS_XV[i,j,k,Bx_])^2 +
            (LS_XV[i,j,k,By_] + RS_XV[i,j,k,By_])^2 +
            (LS_XV[i,j,k,Bz_] + RS_XV[i,j,k,Bz_])^2 ) /
            (2.0*(LS_XV[i,j,k,Rho_] + RS_XV[i,j,k,Rho_]))
         Can2_XF = ( (LS_XV[i,j,k,Bx_] .+ RS_XV[i,j,k,Bx_])^2 ) /
            (2.0 *(LS_XV[i,j,k,Rho_] + RS_XV[i,j,k,Rho_]))
         Cmax_XF[i,j,k] = 0.5 * abs(LS_XV[i,j,k,Ux_]/LS_XV[i,j,k,Rho_] +
            RS_XV[i,j,k,Ux_]/RS_XV[i,j,k,Rho_]) + sqrt( 0.5*(Cs2_XF + Ca2_XF +
            sqrt((Cs2_XF + Ca2_XF)^2 - 4*Cs2_XF*Can2_XF)) )
      end

      for k = 1:nK, j = 1:nJ+1, i = 1:nI
         Cs2_YF = γ*(LS_YV[i,j,k,P_] + RS_YV[i,j,k,P_]) /
            (LS_YV[i,j,k,Rho_] + RS_YV[i,j,k,Rho_])
         Ca2_YF = ( (LS_YV[i,j,k,Bx_] + RS_YV[i,j,k,Bx_])^2 +
            (LS_YV[i,j,k,By_] + RS_YV[i,j,k,By_])^2 +
            (LS_YV[i,j,k,Bz_] + RS_YV[i,j,k,Bz_])^2 ) /
            (2.0 *(LS_YV[i,j,k,Rho_] + RS_YV[i,j,k,Rho_]))
         Can2_YF = ( (LS_YV[i,j,k,By_] + RS_YV[i,j,k,By_])^2 ) /
            (2.0 *(LS_YV[i,j,k,Rho_] + RS_YV[i,j,k,Rho_]))
         Cmax_YF[i,j,k] = 0.5 * abs(LS_YV[i,j,k,Uy_]/LS_YV[i,j,k,Rho_] +
            RS_YV[i,j,k,Uy_]/RS_YV[i,j,k,Rho_]) + sqrt( 0.5*(Cs2_YF + Ca2_YF +
            sqrt((Cs2_YF + Ca2_YF)^2 - 4*Cs2_YF*Can2_YF)) )
      end

      for k = 1:nK+1, j = 1:nJ, i = 1:nI
         Cs2_ZF = γ*(LS_ZV[i,j,k,P_] + RS_ZV[i,j,k,P_]) /
            (LS_ZV[i,j,k,Rho_] + RS_ZV[i,j,k,Rho_])
         Ca2_ZF = ( (LS_ZV[i,j,k,Bx_] + RS_ZV[i,j,k,Bx_])^2 +
            (LS_ZV[i,j,k,By_] + RS_ZV[i,j,k,By_])^2 +
            (LS_ZV[i,j,k,Bz_] + RS_ZV[i,j,k,Bz_])^2 ) /
            (2.0 *(LS_ZV[i,j,k,Rho_] + RS_ZV[i,j,k,Rho_]))
         Can2_ZF = ( (LS_ZV[i,j,k,Bz_] + RS_ZV[i,j,k,Bz_])^2 ) /
            (2.0 *(LS_ZV[i,j,k,Rho_] + RS_ZV[i,j,k,Rho_]))
         Cmax_ZF[i,j,k] = 0.5 * abs(LS_ZV[i,j,k,Uz_]/LS_ZV[i,j,k,Rho_] +
            RS_ZV[i,j,k,Uz_]/RS_ZV[i,j,k,Rho_]) + sqrt( 0.5*(Cs2_ZF + Ca2_ZF +
            sqrt((Cs2_ZF + Ca2_ZF)^2 - 4*Cs2_ZF*Can2_ZF)) )
      end

   end

   function get_speed_maxmin(faceValue::FaceState)

   end

   end
