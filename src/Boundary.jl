module Boundary

export set_cell_boundary!

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_

function set_cell_boundary!(param::Param, state_GV)

   nG, TypeBc = param.nG, param.TypeBc
   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax
   iMinAll, iMaxAll = param.iMinAll, param.iMaxAll
   jMinAll, jMaxAll = param.jMinAll, param.jMaxAll
   kMinAll, kMaxAll = param.kMinAll, param.kMaxAll
   nVar = param.nVar

   for (iBc,TypeBc) in enumerate(TypeBc)
      if TypeBc == "periodic"
         if iBc == 1
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=jMin:jMax, i=1:nG
               state_GV[i,j,k,iVar] = state_GV[iMax-nG+i,j,k,iVar]
            end
         elseif iBc == 2
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=jMin:jMax, i=iMax+1:iMaxAll
               state_GV[i,j,k,iVar] = state_GV[i-iMax,j,k,iVar]
            end
         elseif iBc == 3
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=1:nG, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,jMax-nG+j,k,iVar]
            end
         elseif iBc == 4
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=iMax+1:iMaxAll, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,j-jMax,k,iVar]
            end
         elseif iBc == 5
            @inbounds for iVar=1:nVar, k=1:nG, j=jMin:jMax, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,j,kMax-nG+k,iVar]
            end
         elseif iBc == 6
            @inbounds for iVar=1:nVar, k=kMax+1:kMaxAll, j=jMin:jMax, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,j,k-kMax,iVar]
            end
         end
      elseif TypeBc == "fixed"
         #@show CellState_GV
         # this needs to be changed
         if iBc == 1
            state_GV[1:nG,jMin:jMax,kMin:kMax,:] .= CellState_GV
         elseif iBc == 2
            state_GV[iMax+1:iMaxAll,jMin:jMax,kMin:kMax,:] .= CellState_GV
         elseif iBc == 3
            state_GV[iMin:iMax,1:nG,kMin:kMax,:] .= CellState_GV
         elseif iBc == 4
            state_GV[iMin:iMax,jMax+1:jMaxAll,kMin:kMax,:] .= CellState_GV
         elseif iBc == 5
            state_GV[iMin:iMax,jMin:jMax,1:nG,:] .= CellState_GV
         elseif iBc == 6
            state_GV[iMin:iMax,jMin:jMax,kMax+1:kMaxAll,:] .= CellState_GV
         end
      elseif TypeBc == "float"
         if iBc == 1
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=jMin:jMax, i=1:nG
               state_GV[i,j,k,iVar] = state_GV[iMin,j,k,iVar]
            end
         elseif iBc == 2
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=jMin:jMax, i=iMax+1:iMaxAll
               state_GV[i,j,k,iVar] = state_GV[iMax,j,k,iVar]
            end
         elseif iBc == 3
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=1:nG, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,jMin,k,iVar]
            end
         elseif iBc == 4
            @inbounds for iVar=1:nVar, k=kMin:kMax, j=jMax+1:jMaxAll, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,jMax,k,iVar]
            end
         elseif iBc == 5
            @inbounds for iVar=1:nVar, k=1:nG, j=jMin:jMax, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,j,kMin,iVar]
            end
         elseif iBc == 6
            @inbounds for iVar=1:nVar, k=kMax+1:kMaxAll, j=jMin:jMax, i=iMin:iMax
               state_GV[i,j,k,iVar] = state_GV[i,j,kMax,iVar]
            end
         end
      end
   end

   return nothing
end

end
