# Test of view usage

const Rho_       = 1
const Ux_        = 2
const Uy_        = 3
const Uz_        = 4
const Bx_        = 5
const By_        = 6
const Bz_        = 7
const P_         = 8

function flux_test()

   nI, nJ, nK, nVar = 40, 40, 40, 8
   nG = 1

   state_GV = zeros(nI+2*nG,nJ+2*nG,nK+2*nG,nVar)
   Flux_XV = zeros(nI+1,nJ,nK,nVar)

   iMin,iMax = 1+nG, nI+nG
   jMin,jMax = 1+nG, nJ+nG
   kMin,kMax = 1+nG, nK+nG

   LS_XV = @view state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
   RS_XV = @view state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:]

   @time for k = 1:nK, j = 1:nJ, i = 1:nI+1
      bL = LS_XV[i,j,k,Bx_]^2 + LS_XV[i,j,k,By_]^2 + LS_XV[i,j,k,Bz_]^2
      Flux_XV[i,j,k,Ux_] = LS_XV[i,j,k,Ux_]^2 / LS_XV[i,j,k,Rho_] +
         LS_XV[i,j,k,P_] + 0.5*bL - LS_XV[i,j,k,Bx_]^2
   end


   @time @. Flux_XV[:,:,:,Ux_] = LS_XV[:,:,:,Ux_]^2 / LS_XV[:,:,:,Rho_] + LS_XV[:,:,:,P_] +
      0.5*(LS_XV[:,:,:,Bx_]^2 + LS_XV[:,:,:,By_]^2 + LS_XV[:,:,:,Bz_]^2) -
      LS_XV[:,:,:,Bx_]^2

   LS_XV = state_GV[iMin-1:iMax,jMin:jMax,kMin:kMax,:]
   RS_XV = state_GV[iMin:iMax+1,jMin:jMax,kMin:kMax,:]


   @time for k = 1:nK, j = 1:nJ, i = 1:nI+1
      bL = LS_XV[i,j,k,Bx_]^2 + LS_XV[i,j,k,By_]^2 + LS_XV[i,j,k,Bz_]^2
      Flux_XV[i,j,k,Ux_] = LS_XV[i,j,k,Ux_]^2 / LS_XV[i,j,k,Rho_] +
         LS_XV[i,j,k,P_] + 0.5*bL - LS_XV[i,j,k,Bx_]^2
   end

   @time @. Flux_XV[:,:,:,Ux_] = LS_XV[:,:,:,Ux_]^2 / LS_XV[:,:,:,Rho_] + LS_XV[:,:,:,P_] +
      0.5*(LS_XV[:,:,:,Bx_]^2 + LS_XV[:,:,:,By_]^2 + LS_XV[:,:,:,Bz_]^2) -
      LS_XV[:,:,:,Bx_]^2

   return
end
