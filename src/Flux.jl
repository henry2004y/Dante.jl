module Flux

export calc_face_flux!, init_flux, SpeedFlux, FaceFlux, FaceFluxLR

using ..Parameters: Param, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, E_, U_, B_
using ..Parameters: γ
using ..FaceValue: FaceState

abstract type SpeedFlux end

struct SpeedFluxMax <: SpeedFlux
   Cmax_XF::Array{Float64,3}
   Cmax_YF::Array{Float64,3}
   Cmax_ZF::Array{Float64,3}
end

struct SpeedFluxMinMax <: SpeedFlux
   Cmax_XF::Array{Float64,3}
   Cmax_YF::Array{Float64,3}
   Cmax_ZF::Array{Float64,3}
   Cmin_XF::Array{Float64,3}
   Cmin_YF::Array{Float64,3}
   Cmin_ZF::Array{Float64,3}
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

   faceFluxLR = FaceFluxLR(LFlux_XV,RFlux_XV,LFlux_YV,
      RFlux_YV,LFlux_ZV,RFlux_ZV)
   faceFlux = FaceFlux(Flux_XV, Flux_YV, Flux_ZV)
   if param.Scheme == "Rusanov"
      speedFlux = SpeedFluxMax(Cmax_XF, Cmax_YF, Cmax_ZF)
   elseif param.Scheme == "HLLE"
      Cmin_XF = similar(Cmax_XF)
      Cmin_YF = similar(Cmax_YF)
      Cmin_ZF = similar(Cmax_ZF)
      speedFlux = SpeedFluxMinMax(Cmax_XF,Cmax_YF,Cmax_ZF,
         Cmin_XF,Cmin_YF,Cmin_ZF)
   end

   return faceFluxLR, faceFlux, speedFlux
end

function calc_face_flux!(param::Param, faceValue::FaceState,
   faceFlux::FaceFlux, speedFlux::SpeedFlux, faceFluxLR::FaceFluxLR)

   get_physical_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

   add_numerical_flux!(param, faceValue, faceFlux, speedFlux, faceFluxLR)

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
   @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
      LFlux_XV[i,j,k,Rho_] = LState_XV[i,j,k,Ux_]
      RFlux_XV[i,j,k,Rho_] = RState_XV[i,j,k,Ux_]
   end

   @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
      LFlux_YV[i,j,k,Rho_] = LState_YV[i,j,k,Uy_]
      RFlux_YV[i,j,k,Rho_] = RState_YV[i,j,k,Uy_]
   end

   @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
      LFlux_ZV[i,j,k,Rho_] = LState_ZV[i,j,k,Uz_]
      RFlux_ZV[i,j,k,Rho_] = RState_ZV[i,j,k,Uz_]
   end

   # Momentum flux
   @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
      bL = LState_XV[i,j,k,Bx_]^2 + LState_XV[i,j,k,By_]^2 +
         LState_XV[i,j,k,Bz_]^2
      bR = RState_XV[i,j,k,Bx_]^2 + RState_XV[i,j,k,By_]^2 +
         RState_XV[i,j,k,Bz_]^2

      LFlux_XV[i,j,k,Ux_] = LState_XV[i,j,k,Ux_]^2 / LState_XV[i,j,k,Rho_] +
         LState_XV[i,j,k,P_] + 0.5*bL - LState_XV[i,j,k,Bx_]^2
      RFlux_XV[i,j,k,Ux_] = RState_XV[i,j,k,Ux_]^2 / RState_XV[i,j,k,Rho_] +
         RState_XV[i,j,k,P_] + 0.5*bR - RState_XV[i,j,k,Bx_]^2
      LFlux_XV[i,j,k,Uy_] = LState_XV[i,j,k,Ux_] * LState_XV[i,j,k,Uy_] /
         LState_XV[i,j,k,Rho_] - LState_XV[i,j,k,Bx_]*LState_XV[i,j,k,By_]
      RFlux_XV[i,j,k,Uy_] = RState_XV[i,j,k,Ux_] * RState_XV[i,j,k,Uy_] /
         RState_XV[i,j,k,Rho_] - RState_XV[i,j,k,Bx_]*RState_XV[i,j,k,By_]
      LFlux_XV[i,j,k,Uz_] = LState_XV[i,j,k,Ux_] * LState_XV[i,j,k,Uz_] /
         LState_XV[i,j,k,Rho_] - LState_XV[i,j,k,Bx_]*LState_XV[i,j,k,Bz_]
      RFlux_XV[i,j,k,Uz_] = RState_XV[i,j,k,Ux_] * RState_XV[i,j,k,Uz_] /
         RState_XV[i,j,k,Rho_] - RState_XV[i,j,k,Bx_]*RState_XV[i,j,k,Bz_]
   end

   @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
      bL = LState_YV[i,j,k,Bx_]^2 + LState_YV[i,j,k,By_]^2 +
         LState_YV[i,j,k,Bz_]^2
      bR = RState_YV[i,j,k,Bx_]^2 + RState_YV[i,j,k,By_]^2 +
         RState_YV[i,j,k,Bz_]^2

      LFlux_YV[i,j,k,Ux_] = LState_YV[i,j,k,Uy_] * LState_YV[i,j,k,Ux_] /
         LState_YV[i,j,k,Rho_] - LState_YV[i,j,k,By_]*LState_YV[i,j,k,Bx_]
      RFlux_YV[i,j,k,Ux_] = RState_YV[i,j,k,Uy_] * RState_YV[i,j,k,Ux_] /
         RState_YV[i,j,k,Rho_] - RState_YV[i,j,k,By_]*RState_YV[i,j,k,Bx_]
      LFlux_YV[i,j,k,Uy_] = LState_YV[i,j,k,Uy_]^2 / LState_YV[i,j,k,Rho_] +
         LState_YV[i,j,k,P_] + 0.5*bL - LState_YV[i,j,k,By_]^2
      RFlux_YV[i,j,k,Uy_] = RState_YV[i,j,k,Uy_]^2 / RState_YV[i,j,k,Rho_] +
         RState_YV[i,j,k,P_] + 0.5*bR - RState_YV[i,j,k,By_]^2
      LFlux_YV[i,j,k,Uz_] = LState_YV[i,j,k,Uy_] * LState_YV[i,j,k,Uz_] /
         LState_YV[i,j,k,Rho_] - LState_YV[i,j,k,Bx_]*LState_YV[i,j,k,Bz_]
      RFlux_YV[i,j,k,Uz_] = RState_YV[i,j,k,Uy_] * RState_YV[i,j,k,Uz_] /
         RState_YV[i,j,k,Rho_] - RState_YV[i,j,k,Bx_]*RState_YV[i,j,k,Bz_]
   end

   @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
      bL = LState_ZV[i,j,k,Bx_]^2 + LState_ZV[i,j,k,By_]^2 +
         LState_ZV[i,j,k,Bz_]^2
      bR = RState_ZV[i,j,k,Bx_]^2 + RState_ZV[i,j,k,By_]^2 +
         RState_ZV[i,j,k,Bz_]^2

      LFlux_ZV[i,j,k,Ux_] = LState_ZV[i,j,k,Uz_] * LState_ZV[i,j,k,Ux_] /
         LState_ZV[i,j,k,Rho_] - LState_ZV[i,j,k,Bz_]*LState_ZV[i,j,k,Bx_]
      RFlux_ZV[i,j,k,Ux_] = RState_ZV[i,j,k,Uz_] * RState_ZV[i,j,k,Ux_] /
         RState_ZV[i,j,k,Rho_] - RState_ZV[i,j,k,Bz_]*RState_ZV[i,j,k,Bx_]
      LFlux_ZV[i,j,k,Uy_] = LState_ZV[i,j,k,Uz_] * LState_ZV[i,j,k,Uy_] /
         LState_ZV[i,j,k,Rho_] - LState_ZV[i,j,k,Bz_]*LState_ZV[i,j,k,By_]
      RFlux_ZV[i,j,k,Uy_] = RState_ZV[i,j,k,Uz_] * RState_ZV[i,j,k,Uy_] /
         RState_ZV[i,j,k,Rho_] - RState_ZV[i,j,k,Bz_]*RState_ZV[i,j,k,By_]
      LFlux_ZV[i,j,k,Uz_] = LState_ZV[i,j,k,Uz_]^2 / LState_ZV[i,j,k,Rho_] +
         LState_ZV[i,j,k,P_] + 0.5*bL - LState_ZV[i,j,k,Bz_]^2
      RFlux_ZV[i,j,k,Uz_] = RState_ZV[i,j,k,Uz_]^2 / RState_ZV[i,j,k,Rho_] +
         RState_ZV[i,j,k,P_] + 0.5*bR - RState_ZV[i,j,k,Bz_]^2
   end

   # Magnetic flux
   @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
      LFlux_XV[i,j,k,Bx_] = 0.0
      RFlux_XV[i,j,k,Bx_] = 0.0
      LFlux_XV[i,j,k,By_] = LState_XV[i,j,k,Ux_]*LState_XV[i,j,k,By_] -
         LState_XV[i,j,k,Bx_]*LState_XV[i,j,k,Uy_]
      RFlux_XV[i,j,k,By_] = RState_XV[i,j,k,Ux_]*RState_XV[i,j,k,By_] -
         RState_XV[i,j,k,Bx_]*RState_XV[i,j,k,Uy_]
      LFlux_XV[i,j,k,Bz_] = LState_XV[i,j,k,Ux_]*LState_XV[i,j,k,Bz_] -
         LState_XV[i,j,k,Bx_]*LState_XV[i,j,k,Uz_]
      RFlux_XV[i,j,k,Bz_] = RState_XV[i,j,k,Ux_]*RState_XV[i,j,k,Bz_] -
         RState_XV[i,j,k,Bx_]*RState_XV[i,j,k,Uz_]
   end

   @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
      LFlux_YV[i,j,k,Bx_] = LState_YV[i,j,k,Uy_]*LState_YV[i,j,k,Bx_] -
         LState_YV[i,j,k,By_]*LState_YV[i,j,k,Ux_]
      RFlux_YV[i,j,k,Bx_] = RState_YV[i,j,k,Uy_]*RState_YV[i,j,k,Bx_] -
         RState_YV[i,j,k,By_]*RState_YV[i,j,k,Ux_]
      LFlux_YV[i,j,k,By_] = 0.0
      RFlux_YV[i,j,k,By_] = 0.0
      LFlux_YV[i,j,k,Bz_] = LState_YV[i,j,k,Uy_]*LState_YV[i,j,k,Bz_] -
         LState_YV[i,j,k,By_]*LState_YV[i,j,k,Uz_]
      RFlux_YV[i,j,k,Bz_] = RState_YV[i,j,k,Uy_]*RState_YV[i,j,k,Bz_] -
         RState_YV[i,j,k,By_]*RState_YV[i,j,k,Uz_]
   end

   @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
      LFlux_ZV[i,j,k,Bx_] = LState_ZV[i,j,k,Uz_]*LState_ZV[i,j,k,Bx_] -
         LState_ZV[i,j,k,Bz_]*LState_ZV[i,j,k,Ux_]
      RFlux_ZV[i,j,k,Bx_] = RState_ZV[i,j,k,Uz_]*RState_ZV[i,j,k,Bx_] -
         RState_ZV[i,j,k,Bz_]*RState_ZV[i,j,k,Ux_]
      LFlux_ZV[i,j,k,By_] = LState_ZV[i,j,k,Uz_]*LState_ZV[i,j,k,By_] -
         LState_ZV[i,j,k,Bz_]*LState_ZV[i,j,k,Uy_]
      RFlux_ZV[i,j,k,By_] = RState_ZV[i,j,k,Uz_]*RState_ZV[i,j,k,By_] -
         RState_ZV[i,j,k,Bz_]*RState_ZV[i,j,k,Uy_]
      LFlux_ZV[i,j,k,Bz_] = 0.0
      RFlux_ZV[i,j,k,Bz_] = 0.0
   end

   # Pressure flux / energy flux
   if !param.UseConservative
      @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
         LFlux_XV[i,j,k,P_] = LState_XV[i,j,k,P_]*LState_XV[i,j,k,Ux_]/
            LState_XV[i,j,k,Rho_]
         RFlux_XV[i,j,k,P_] = RState_XV[i,j,k,P_]*RState_XV[i,j,k,Ux_]/
            RState_XV[i,j,k,Rho_]
      end

      @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
         LFlux_YV[i,j,k,P_] = LState_YV[i,j,k,P_]*LState_YV[i,j,k,Uy_]/
            LState_YV[i,j,k,Rho_]
         RFlux_YV[i,j,k,P_] = RState_YV[i,j,k,P_]*RState_YV[i,j,k,Uy_]/
            RState_YV[i,j,k,Rho_]
      end

      @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
         LFlux_ZV[i,j,k,P_] = LState_ZV[i,j,k,P_]*LState_ZV[i,j,k,Uz_]/
            LState_ZV[i,j,k,Rho_]
         RFlux_ZV[i,j,k,P_] = RState_ZV[i,j,k,P_]*RState_ZV[i,j,k,Uz_]/
            RState_ZV[i,j,k,Rho_]
      end
   else
      # Currently use the same index for pressure/energy
      @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
         u = LState_XV[i,j,k,Ux_]^2 + LState_XV[i,j,k,Uy_]^2 +
            LState_XV[i,j,k,Uz_]^2
         b = LState_XV[i,j,k,Bx_]^2 + LState_XV[i,j,k,By_]^2 +
            LState_XV[i,j,k,Bz_]^2
         ub= LState_XV[i,j,k,Ux_]*LState_XV[i,j,k,Bx_] +
            LState_XV[i,j,k,Uy_]*LState_XV[i,j,k,By_] +
            LState_XV[i,j,k,Uz_]*LState_XV[i,j,k,Bz_]
         LFlux_XV[i,j,k,E_] = LState_XV[i,j,k,Ux_]/LState_XV[i,j,k,Rho_]*
            ((LState_XV[i,j,k,P_]/(γ-1.0) + 0.5/LState_XV[i,j,k,Rho_]*u +
            0.5*b + LState_XV[i,j,k,P_] + 0.5*b)) - ub *
            LState_XV[i,j,k,Bx_]
      end

      @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
         u = RState_XV[i,j,k,Ux_]^2 + RState_XV[i,j,k,Uy_]^2 +
            RState_XV[i,j,k,Uz_]^2
         b = RState_XV[i,j,k,Bx_]^2 + RState_XV[i,j,k,By_]^2 +
            RState_XV[i,j,k,Bz_]^2
         ub= RState_XV[i,j,k,Ux_]*RState_XV[i,j,k,Bx_] +
            RState_XV[i,j,k,Uy_]*RState_XV[i,j,k,By_] +
            RState_XV[i,j,k,Uz_]*RState_XV[i,j,k,Bz_]
         RFlux_XV[i,j,k,E_] = RState_XV[i,j,k,Ux_]/RState_XV[i,j,k,Rho_]*
            ((RState_XV[i,j,k,P_]/(γ-1.0) + 0.5/RState_XV[i,j,k,Rho_]*u +
            0.5*b + RState_XV[i,j,k,P_] + 0.5*b)) - ub *
            RState_XV[i,j,k,Bx_]
      end

      @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
         u = LState_YV[i,j,k,Ux_]^2 + LState_YV[i,j,k,Uy_]^2 +
            LState_YV[i,j,k,Uz_]^2
         b = LState_YV[i,j,k,Bx_]^2 + LState_YV[i,j,k,By_]^2 +
            LState_YV[i,j,k,Bz_]^2
         ub= LState_YV[i,j,k,Ux_]*LState_YV[i,j,k,Bx_] +
            LState_YV[i,j,k,Uy_]*LState_YV[i,j,k,By_] +
            LState_YV[i,j,k,Uz_]*LState_YV[i,j,k,Bz_]
         LFlux_YV[i,j,k,E_] = LState_YV[i,j,k,Uy_]/LState_YV[i,j,k,Rho_]*
            ((LState_YV[i,j,k,P_]/(γ-1.0) + 0.5/LState_YV[i,j,k,Rho_]*u +
            0.5*b + LState_YV[i,j,k,P_] + 0.5*b)) - ub *
            LState_YV[i,j,k,Bx_]
      end

      @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
         u = RState_YV[i,j,k,Ux_]^2 + RState_YV[i,j,k,Uy_]^2 +
            RState_YV[i,j,k,Uz_]^2
         b = RState_YV[i,j,k,Bx_]^2 + RState_YV[i,j,k,By_]^2 +
            RState_YV[i,j,k,Bz_]^2
         ub= RState_YV[i,j,k,Ux_]*RState_YV[i,j,k,Bx_] +
            RState_YV[i,j,k,Uy_]*RState_YV[i,j,k,By_] +
            RState_YV[i,j,k,Uz_]*RState_YV[i,j,k,Bz_]
         RFlux_YV[i,j,k,E_] = RState_YV[i,j,k,Uy_]/RState_YV[i,j,k,Rho_]*
            ((RState_YV[i,j,k,P_]/(γ-1.0) + 0.5/RState_YV[i,j,k,Rho_]*u +
            0.5*b + RState_YV[i,j,k,P_] + 0.5*b)) - ub *
            RState_YV[i,j,k,Bx_]
      end

      @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
         u = LState_ZV[i,j,k,Ux_]^2 + LState_ZV[i,j,k,Uy_]^2 +
            LState_ZV[i,j,k,Uz_]^2
         b = LState_ZV[i,j,k,Bx_]^2 + LState_ZV[i,j,k,By_]^2 +
            LState_ZV[i,j,k,Bz_]^2
         ub= LState_ZV[i,j,k,Ux_]*LState_ZV[i,j,k,Bx_] +
            LState_ZV[i,j,k,Uy_]*LState_ZV[i,j,k,By_] +
            LState_ZV[i,j,k,Uz_]*LState_ZV[i,j,k,Bz_]
         LFlux_ZV[i,j,k,E_] = LState_ZV[i,j,k,Uz_]/LState_ZV[i,j,k,Rho_]*
            ((LState_ZV[i,j,k,P_]/(γ-1.0) + 0.5/LState_ZV[i,j,k,Rho_]*u +
            0.5*b + LState_ZV[i,j,k,P_] + 0.5*b)) - ub *
            LState_ZV[i,j,k,Bx_]
      end

      @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
         u = RState_ZV[i,j,k,Ux_]^2 + RState_ZV[i,j,k,Uy_]^2 +
            RState_ZV[i,j,k,Uz_]^2
         b = RState_ZV[i,j,k,Bx_]^2 + RState_ZV[i,j,k,By_]^2 +
            RState_ZV[i,j,k,Bz_]^2
         ub= RState_ZV[i,j,k,Ux_]*RState_ZV[i,j,k,Bx_] +
            RState_ZV[i,j,k,Uy_]*RState_ZV[i,j,k,By_] +
            RState_ZV[i,j,k,Uz_]*RState_ZV[i,j,k,Bz_]
         RFlux_ZV[i,j,k,E_] = RState_ZV[i,j,k,Uz_]/RState_ZV[i,j,k,Rho_]*
            ((RState_ZV[i,j,k,P_]/(γ-1.0) + 0.5/RState_ZV[i,j,k,Rho_]*u +
            0.5*b + RState_ZV[i,j,k,P_] + 0.5*b)) - ub *
            RState_ZV[i,j,k,Bx_]
      end
   end

   # Collect all the physical fluxes
   @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
      Flux_XV[i,j,k,iVar] = 0.5*(LFlux_XV[i,j,k,iVar] + RFlux_XV[i,j,k,iVar])
   end
   @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
      Flux_YV[i,j,k,iVar] = 0.5*(LFlux_YV[i,j,k,iVar] + RFlux_YV[i,j,k,iVar])
   end
   @inbounds for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
      Flux_ZV[i,j,k,iVar] = 0.5*(LFlux_ZV[i,j,k,iVar] + RFlux_ZV[i,j,k,iVar])
   end

end

function add_numerical_flux!(param::Param, faceValue::FaceState,
   faceFlux::FaceFlux, speedFlux::SpeedFlux, faceFluxLR::FaceFluxLR)

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
         @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
            Flux_XV[i,j,k,iVar] -= 0.5*Cmax_XF[i,j,k]*(RState_XV[i,j,k,iVar] -
               LState_XV[i,j,k,iVar])
         end

         @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
            Flux_YV[i,j,k,iVar] -= 0.5*Cmax_YF[i,j,k]*(RState_YV[i,j,k,iVar] -
               LState_YV[i,j,k,iVar])
         end

         @inbounds for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
            Flux_ZV[i,j,k,iVar] -= 0.5*Cmax_ZF[i,j,k]*(RState_ZV[i,j,k,iVar] -
               LState_ZV[i,j,k,iVar])
         end
      else
         # If I solve energy equation instead of pressure, there's
         # duplicate calculation above, even though the expression
         # looks compact. That's why I use an if-else statement
         @inbounds for iVar = Rho_:Bz_, k = 1:nK, j = 1:nJ, i = 1:nI+1
            Flux_XV[i,j,k,iVar] -= 0.5*Cmax_XF[i,j,k]*(RState_XV[i,j,k,iVar] -
               LState_XV[i,j,k,iVar])
         end

         @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
            uL = LState_XV[i,j,k,Ux_]^2 + LState_XV[i,j,k,Uy_]^2 +
               LState_XV[i,j,k,Uz_]^2
            uR = RState_XV[i,j,k,Ux_]^2 + RState_XV[i,j,k,Uy_]^2 +
               RState_XV[i,j,k,Uz_]^2
            bL = LState_XV[i,j,k,Bx_]^2 + LState_XV[i,j,k,By_]^2 +
               LState_XV[i,j,k,Bz_]^2
            bR = RState_XV[i,j,k,Bx_]^2 + RState_XV[i,j,k,By_]^2 +
               RState_XV[i,j,k,Bz_]^2

            Flux_XV[i,j,k,E_] -= 0.5*Cmax_XF[i,j,k]*(
               (RState_XV[i,j,k,P_]/(γ-1) + 0.5/RState_XV[i,j,k,Rho_]*uR
               + 0.5*bR) -
               (LState_XV[i,j,k,P_]/(γ-1) + 0.5/LState_XV[i,j,k,Rho_]*uL
               + 0.5*bL))
         end

         @inbounds for iVar = Rho_:Bz_, k = 1:nK, j = 1:nJ+1, i = 1:nI
            Flux_YV[i,j,k,iVar] -= 0.5*Cmax_YF[i,j,k]*(RState_YV[i,j,k,iVar] -
               LState_YV[i,j,k,iVar])
         end

         @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
            uL = LState_YV[i,j,k,Ux_]^2 + LState_YV[i,j,k,Uy_]^2 +
               LState_YV[i,j,k,Uz_]^2
            uR = RState_YV[i,j,k,Ux_]^2 + RState_YV[i,j,k,Uy_]^2 +
               RState_YV[i,j,k,Uz_]^2
            bL = LState_YV[i,j,k,Bx_]^2 + LState_YV[i,j,k,By_]^2 +
               LState_YV[i,j,k,Bz_]^2
            bR = RState_YV[i,j,k,Bx_]^2 + RState_YV[i,j,k,By_]^2 +
               RState_YV[i,j,k,Bz_]^2

            Flux_YV[i,j,k,E_] -= 0.5*Cmax_YF[i,j,k]*(
               (RState_YV[i,j,k,P_]/(γ-1) + 0.5/RState_YV[i,j,k,Rho_]*uR
               + 0.5*bR) -
               (LState_YV[i,j,k,P_]/(γ-1) + 0.5/LState_YV[i,j,k,Rho_]*uL
               + 0.5*bL))
         end

         @inbounds for iVar = Rho_:Bz_, k = 1:nK+1, j = 1:nJ, i = 1:nI
            Flux_ZV[i,j,k,iVar] -= 0.5*Cmax_ZF[i,j,k]*(RState_ZV[i,j,k,iVar] -
               LState_ZV[i,j,k,iVar])
         end

         @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
            uL = LState_ZV[i,j,k,Ux_]^2 + LState_ZV[i,j,k,Uy_]^2 +
               LState_ZV[i,j,k,Uz_]^2
            uR = RState_ZV[i,j,k,Ux_]^2 + RState_ZV[i,j,k,Uy_]^2 +
               RState_ZV[i,j,k,Uz_]^2
            bL = LState_ZV[i,j,k,Bx_]^2 + LState_ZV[i,j,k,By_]^2 +
               LState_ZV[i,j,k,Bz_]^2
            bR = RState_ZV[i,j,k,Bx_]^2 + RState_ZV[i,j,k,By_]^2 +
               RState_ZV[i,j,k,Bz_]^2

            Flux_ZV[i,j,k,E_] -= 0.5*Cmax_ZF[i,j,k]*(
               (RState_ZV[i,j,k,P_]/(γ-1) + 0.5/RState_ZV[i,j,k,Rho_]*uR
               + 0.5*bR) -
               (LState_ZV[i,j,k,P_]/(γ-1) + 0.5/LState_ZV[i,j,k,Rho_]*uL
               + 0.5*bL))
         end

      end
   elseif param.Scheme == "HLLE"
      Cmin_XF = speedFlux.Cmin_XF
      Cmin_YF = speedFlux.Cmin_YF
      Cmin_ZF = speedFlux.Cmin_ZF

      LFlux_XV = faceFluxLR.LFlux_XV
      RFlux_XV = faceFluxLR.RFlux_XV
      LFlux_YV = faceFluxLR.LFlux_YV
      RFlux_YV = faceFluxLR.RFlux_YV
      LFlux_ZV = faceFluxLR.LFlux_ZV
      RFlux_ZV = faceFluxLR.RFlux_ZV

      get_speed_maxmin!(param, faceValue, speedFlux)

      if ~param.UseConservative
         @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ, i = 1:nI+1
            Flux_XV[i,j,k,iVar] -= 0.5*(Cmax_XF[i,j,k] + Cmin_XF[i,j,k])/
               (Cmax_XF[i,j,k] - Cmin_XF[i,j,k])*
               (RFlux_XV[i,j,k,iVar] - LFlux_XV[i,j,k,iVar]) -
               Cmax_XF[i,j,k]*Cmin_XF[i,j,k]/(Cmax_XF[i,j,k] - Cmin_XF[i,j,k])*
               (RState_XV[i,j,k,iVar] - LState_XV[i,j,k,iVar])
         end

         @inbounds for iVar = 1:nVar, k = 1:nK, j = 1:nJ+1, i = 1:nI
            Flux_YV[i,j,k,iVar] -= 0.5*(Cmax_YF[i,j,k] + Cmin_YF[i,j,k])/
               (Cmax_YF[i,j,k] - Cmin_YF[i,j,k])*
               (RFlux_YV[i,j,k,iVar] - LFlux_YV[i,j,k,iVar]) -
               Cmax_YF[i,j,k]*Cmin_YF[i,j,k]/(Cmax_YF[i,j,k] - Cmin_YF[i,j,k])*
               (RState_YV[i,j,k,iVar] - LState_YV[i,j,k,iVar])
         end

         @inbounds for iVar = 1:nVar, k = 1:nK+1, j = 1:nJ, i = 1:nI
            Flux_ZV[i,j,k,iVar] -= 0.5*(Cmax_ZF[i,j,k] + Cmin_ZF[i,j,k])/
               (Cmax_ZF[i,j,k] - Cmin_ZF[i,j,k])*
               (RFlux_ZV[i,j,k,iVar] - LFlux_ZV[i,j,k,iVar]) -
               Cmax_ZF[i,j,k]*Cmin_ZF[i,j,k]/(Cmax_ZF[i,j,k] - Cmin_ZF[i,j,k])*
               (RState_ZV[i,j,k,iVar] - LState_ZV[i,j,k,iVar])
         end
      else
         # If I solve energy equation instead of pressure, there's
         # duplicate calculation above, even though the expression
         # looks compact. That's why I use an if-else statement.

         @inbounds for iVar = Rho_:Bz_, k = 1:nK, j = 1:nJ, i = 1:nI+1
            Flux_XV[i,j,k,iVar] -= 0.5*(Cmax_XF[i,j,k] + Cmin_XF[i,j,k])/
               (Cmax_XF[i,j,k] - Cmin_XF[i,j,k])*
               (RFlux_XV[i,j,k,iVar] - LFlux_XV[i,j,k,iVar]) -
               Cmax_XF[i,j,k]*Cmin_XF[i,j,k]/(Cmax_XF[i,j,k] - Cmin_XF[i,j,k])*
               (RState_XV[i,j,k,iVar] - LState_XV[i,j,k,iVar])
         end

         @inbounds for k = 1:nK, j = 1:nJ, i = 1:nI+1
            uL = LState_XV[i,j,k,Ux_]^2 + LState_XV[i,j,k,Uy_]^2 +
               LState_XV[i,j,k,Uz_]^2
            uR = RState_XV[i,j,k,Ux_]^2 + RState_XV[i,j,k,Uy_]^2 +
               RState_XV[i,j,k,Uz_]^2
            bL = LState_XV[i,j,k,Bx_]^2 + LState_XV[i,j,k,By_]^2 +
               LState_XV[i,j,k,Bz_]^2
            bR = RState_XV[i,j,k,Bx_]^2 + RState_XV[i,j,k,By_]^2 +
               RState_XV[i,j,k,Bz_]^2

            Flux_XV[i,j,k,E_] -= 0.5*(Cmax_XF[i,j,k] + Cmin_XF[i,j,k])/
               (Cmax_XF[i,j,k] - Cmin_XF[i,j,k])*
               (RFlux_XV[i,j,k,E_] - LFlux_XV[i,j,k,E_]) -
               Cmax_XF[i,j,k]*Cmin_XF[i,j,k]/(Cmax_XF[i,j,k] - Cmin_XF[i,j,k])*(
               (RState_XV[i,j,k,P_] / (γ-1) +
               0.5/RState_XV[i,j,k,Rho_]*uR + 0.5*bR) -
               (LState_XV[i,j,k,P_] / (γ-1) +
               0.5/LState_XV[i,j,k,Rho_]*uL + 0.5*bL) )
         end

         @inbounds for iVar = Rho_:Bz_, k = 1:nK, j = 1:nJ+1, i = 1:nI
            Flux_YV[i,j,k,iVar] -= 0.5*(Cmax_YF[i,j,k] + Cmin_YF[i,j,k])/
               (Cmax_YF[i,j,k] - Cmin_YF[i,j,k])*
               (RFlux_YV[i,j,k,iVar] - LFlux_YV[i,j,k,iVar]) -
               Cmax_YF[i,j,k]*Cmin_YF[i,j,k]/(Cmax_YF[i,j,k] - Cmin_YF[i,j,k])*
               (RState_YV[i,j,k,iVar] - LState_YV[i,j,k,iVar])
         end

         @inbounds for k = 1:nK, j = 1:nJ+1, i = 1:nI
            uL = LState_YV[i,j,k,Ux_]^2 + LState_YV[i,j,k,Uy_]^2 +
               LState_YV[i,j,k,Uz_]^2
            uR = RState_YV[i,j,k,Ux_]^2 + RState_YV[i,j,k,Uy_]^2 +
               RState_YV[i,j,k,Uz_]^2
            bL = LState_YV[i,j,k,Bx_]^2 + LState_YV[i,j,k,By_]^2 +
               LState_YV[i,j,k,Bz_]^2
            bR = RState_YV[i,j,k,Bx_]^2 + RState_YV[i,j,k,By_]^2 +
               RState_YV[i,j,k,Bz_]^2

            Flux_YV[i,j,k,E_] -= 0.5*(Cmax_YF[i,j,k] + Cmin_YF[i,j,k])/
               (Cmax_YF[i,j,k] - Cmin_YF[i,j,k])*
               (RFlux_YV[i,j,k,E_] - LFlux_YV[i,j,k,E_]) -
               Cmax_YF[i,j,k]*Cmin_YF[i,j,k]/(Cmax_YF[i,j,k] - Cmin_YF[i,j,k])*(
               (RState_YV[i,j,k,P_] / (γ-1) +
               0.5/RState_YV[i,j,k,Rho_]*uR + 0.5*bR) -
               (LState_YV[i,j,k,P_] / (γ-1) +
               0.5/LState_YV[i,j,k,Rho_]*uL + 0.5*bL) )
         end

         @inbounds for iVar = Rho_:Bz_, k = 1:nK+1, j = 1:nJ, i = 1:nI
            Flux_ZV[i,j,k,iVar] -= 0.5*(Cmax_ZF[i,j,k] + Cmin_ZF[i,j,k])/
               (Cmax_ZF[i,j,k] - Cmin_ZF[i,j,k])*
               (RFlux_ZV[i,j,k,iVar] - LFlux_ZV[i,j,k,iVar]) -
               Cmax_ZF[i,j,k]*Cmin_ZF[i,j,k]/(Cmax_ZF[i,j,k] - Cmin_ZF[i,j,k])*
               (RState_ZV[i,j,k,iVar] - LState_ZV[i,j,k,iVar])
         end

         @inbounds for k = 1:nK+1, j = 1:nJ, i = 1:nI
            uL = LState_ZV[i,j,k,Ux_]^2 + LState_ZV[i,j,k,Uy_]^2 +
               LState_ZV[i,j,k,Uz_]^2
            uR = RState_ZV[i,j,k,Ux_]^2 + RState_ZV[i,j,k,Uy_]^2 +
               RState_ZV[i,j,k,Uz_]^2
            bL = LState_ZV[i,j,k,Bx_]^2 + LState_ZV[i,j,k,By_]^2 +
               LState_ZV[i,j,k,Bz_]^2
            bR = RState_ZV[i,j,k,Bx_]^2 + RState_ZV[i,j,k,By_]^2 +
               RState_ZV[i,j,k,Bz_]^2

            Flux_ZV[i,j,k,E_] -= 0.5*(Cmax_ZF[i,j,k] + Cmin_ZF[i,j,k])/
               (Cmax_ZF[i,j,k] - Cmin_ZF[i,j,k])*
               (RFlux_ZV[i,j,k,E_] - LFlux_ZV[i,j,k,E_]) -
               Cmax_ZF[i,j,k]*Cmin_ZF[i,j,k]/(Cmax_ZF[i,j,k] - Cmin_ZF[i,j,k])*(
               (RState_ZV[i,j,k,P_] / (γ-1) +
               0.5/RState_ZV[i,j,k,Rho_]*uR + 0.5*bR) -
               (LState_ZV[i,j,k,P_] / (γ-1) +
               0.5/LState_ZV[i,j,k,Rho_]*uL + 0.5*bL))
         end
      end

   end
end

"""Calculate the maximum speed in each direction."""
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
      Can2_XF = ( (LS_XV[i,j,k,Bx_] + RS_XV[i,j,k,Bx_])^2 ) /
         (2.0*(LS_XV[i,j,k,Rho_] + RS_XV[i,j,k,Rho_]))
      Cmax_XF[i,j,k] = 0.5 * abs(LS_XV[i,j,k,Ux_]/LS_XV[i,j,k,Rho_] +
         RS_XV[i,j,k,Ux_]/RS_XV[i,j,k,Rho_]) + sqrt( 0.5*(Cs2_XF + Ca2_XF +
         sqrt((Cs2_XF + Ca2_XF)^2 - 4.0*Cs2_XF*Can2_XF)) )
   end

   for k = 1:nK, j = 1:nJ+1, i = 1:nI
      Cs2_YF = γ*(LS_YV[i,j,k,P_] + RS_YV[i,j,k,P_]) /
         (LS_YV[i,j,k,Rho_] + RS_YV[i,j,k,Rho_])
      Ca2_YF = ( (LS_YV[i,j,k,Bx_] + RS_YV[i,j,k,Bx_])^2 +
         (LS_YV[i,j,k,By_] + RS_YV[i,j,k,By_])^2 +
         (LS_YV[i,j,k,Bz_] + RS_YV[i,j,k,Bz_])^2 ) /
         (2.0*(LS_YV[i,j,k,Rho_] + RS_YV[i,j,k,Rho_]))
      Can2_YF = ( (LS_YV[i,j,k,By_] + RS_YV[i,j,k,By_])^2 ) /
         (2.0*(LS_YV[i,j,k,Rho_] + RS_YV[i,j,k,Rho_]))
      Cmax_YF[i,j,k] = 0.5 * abs(LS_YV[i,j,k,Uy_]/LS_YV[i,j,k,Rho_] +
         RS_YV[i,j,k,Uy_]/RS_YV[i,j,k,Rho_]) + sqrt( 0.5*(Cs2_YF + Ca2_YF +
         sqrt((Cs2_YF + Ca2_YF)^2 - 4.0*Cs2_YF*Can2_YF)) )
   end

   for k = 1:nK+1, j = 1:nJ, i = 1:nI
      Cs2_ZF = γ*(LS_ZV[i,j,k,P_] + RS_ZV[i,j,k,P_]) /
         (LS_ZV[i,j,k,Rho_] + RS_ZV[i,j,k,Rho_])
      Ca2_ZF = ( (LS_ZV[i,j,k,Bx_] + RS_ZV[i,j,k,Bx_])^2 +
         (LS_ZV[i,j,k,By_] + RS_ZV[i,j,k,By_])^2 +
         (LS_ZV[i,j,k,Bz_] + RS_ZV[i,j,k,Bz_])^2 ) /
         (2.0*(LS_ZV[i,j,k,Rho_] + RS_ZV[i,j,k,Rho_]))
      Can2_ZF = ( (LS_ZV[i,j,k,Bz_] + RS_ZV[i,j,k,Bz_])^2 ) /
         (2.0*(LS_ZV[i,j,k,Rho_] + RS_ZV[i,j,k,Rho_]))
      Cmax_ZF[i,j,k] = 0.5 * abs(LS_ZV[i,j,k,Uz_]/LS_ZV[i,j,k,Rho_] +
         RS_ZV[i,j,k,Uz_]/RS_ZV[i,j,k,Rho_]) + sqrt( 0.5*(Cs2_ZF + Ca2_ZF +
         sqrt((Cs2_ZF + Ca2_ZF)^2 - 4.0*Cs2_ZF*Can2_ZF)) )
   end

end

function get_speed_maxmin!(param::Param, faceValue::FaceState,
   speedFlux::SpeedFluxMinMax)
   nI, nJ, nK = param.nI, param.nJ, param.nK
   # Aliases
   LS_XV, RS_XV = faceValue.LState_XV, faceValue.RState_XV
   LS_YV, RS_YV = faceValue.LState_YV, faceValue.RState_YV
   LS_ZV, RS_ZV = faceValue.LState_ZV, faceValue.RState_ZV

   Cmax_XF = speedFlux.Cmax_XF
   Cmax_YF = speedFlux.Cmax_YF
   Cmax_ZF = speedFlux.Cmax_ZF
   Cmin_XF = speedFlux.Cmin_XF
   Cmin_YF = speedFlux.Cmin_YF
   Cmin_ZF = speedFlux.Cmin_ZF

   # There must be better ways to do this: LS_XV and RS_XV is just a shift of state_GV,
   # so some repetitive calculation can be avoided!

   for k = 1:nK, j = 1:nJ, i = 1:nI+1
      Cs2_LXF = γ*LS_XV[i,j,k,P_]/LS_XV[i,j,k,Rho_]
      Ca2_LXF = (LS_XV[i,j,k,Bx_]^2 + LS_XV[i,j,k,By_]^2 + LS_XV[i,j,k,Bz_]^2) /
         LS_XV[i,j,k,Rho_]
      Can2_LXF = LS_XV[i,j,k,Bx_]^2 / LS_XV[i,j,k,Rho_]
      u_LXF = LS_XV[i,j,k,Ux_]/LS_XV[i,j,k,Rho_]
      c_LXF = sqrt( 0.5*(Cs2_LXF + Ca2_LXF +
         sqrt((Cs2_LXF + Ca2_LXF)^2 - 4.0*Cs2_LXF*Can2_LXF)) )

      Cs2_RXF = γ*RS_XV[i,j,k,P_]/RS_XV[i,j,k,Rho_]
      Ca2_RXF = (RS_XV[i,j,k,Bx_]^2 + RS_XV[i,j,k,By_]^2 + RS_XV[i,j,k,Bz_]^2) /
         RS_XV[i,j,k,Rho_]
      Can2_RXF = RS_XV[i,j,k,Bx_]^2 / RS_XV[i,j,k,Rho_]
      u_RXF = RS_XV[i,j,k,Ux_]/RS_XV[i,j,k,Rho_]
      c_RXF = sqrt( 0.5*(Cs2_RXF + Ca2_RXF +
         sqrt((Cs2_RXF + Ca2_RXF)^2 - 4.0*Cs2_RXF*Can2_RXF)) )

      Cmax_XF[i,j,k] = max(0, u_LXF+c_LXF, u_RXF+c_RXF)
      Cmin_XF[i,j,k] = min(0, u_LXF-c_LXF, u_RXF-c_RXF)
   end

   for k = 1:nK, j = 1:nJ+1, i = 1:nI
      Cs2_LYF = γ*LS_YV[i,j,k,P_]/LS_YV[i,j,k,Rho_]
      Ca2_LYF = (LS_YV[i,j,k,Bx_]^2 + LS_YV[i,j,k,By_]^2 + LS_YV[i,j,k,Bz_]^2) /
         LS_YV[i,j,k,Rho_]
      Can2_LYF = LS_YV[i,j,k,Bx_]^2 / LS_YV[i,j,k,Rho_]
      u_LYF = LS_YV[i,j,k,Uy_]/LS_YV[i,j,k,Rho_]
      c_LYF = sqrt( 0.5*(Cs2_LYF + Ca2_LYF +
         sqrt((Cs2_LYF + Ca2_LYF)^2 - 4.0*Cs2_LYF*Can2_LYF)) )

      Cs2_RYF = γ*RS_YV[i,j,k,P_]/RS_YV[i,j,k,Rho_]
      Ca2_RYF = (RS_YV[i,j,k,Bx_]^2 + RS_YV[i,j,k,By_]^2 + RS_YV[i,j,k,Bz_]^2) /
         RS_YV[i,j,k,Rho_]
      Can2_RYF = RS_YV[i,j,k,Bx_]^2 / RS_YV[i,j,k,Rho_]
      u_RYF = RS_YV[i,j,k,Ux_]/RS_YV[i,j,k,Rho_]
      c_RYF = sqrt( 0.5*(Cs2_RYF + Ca2_RYF +
         sqrt((Cs2_RYF + Ca2_RYF)^2 - 4.0*Cs2_RYF*Can2_RYF)) )

      Cmax_YF[i,j,k] = max(0, u_LYF+c_LYF, u_RYF+c_RYF)
      Cmin_YF[i,j,k] = min(0, u_LYF-c_LYF, u_RYF-c_RYF)
   end

   for k = 1:nK+1, j = 1:nJ, i = 1:nI
      Cs2_LZF = γ*LS_ZV[i,j,k,P_]/LS_ZV[i,j,k,Rho_]
      Ca2_LZF = (LS_ZV[i,j,k,Bx_]^2 + LS_ZV[i,j,k,By_]^2 + LS_ZV[i,j,k,Bz_]^2) /
         LS_ZV[i,j,k,Rho_]
      Can2_LZF = LS_ZV[i,j,k,Bx_]^2 / LS_ZV[i,j,k,Rho_]
      u_LZF = LS_ZV[i,j,k,Ux_]/LS_ZV[i,j,k,Rho_]
      c_LZF = sqrt( 0.5*(Cs2_LZF + Ca2_LZF +
         sqrt((Cs2_LZF + Ca2_LZF)^2 - 4.0*Cs2_LZF*Can2_LZF)) )

      Cs2_RZF = γ*RS_ZV[i,j,k,P_]/RS_ZV[i,j,k,Rho_]
      Ca2_RZF = (RS_ZV[i,j,k,Bx_]^2 + RS_ZV[i,j,k,By_]^2 + RS_ZV[i,j,k,Bz_]^2) /
         RS_ZV[i,j,k,Rho_]
      Can2_RZF = RS_ZV[i,j,k,Bx_]^2 / RS_ZV[i,j,k,Rho_]
      u_RZF = RS_ZV[i,j,k,Uz_]/RS_ZV[i,j,k,Rho_]
      c_RZF = sqrt( 0.5*(Cs2_RZF + Ca2_RZF +
         sqrt((Cs2_RZF + Ca2_RZF)^2 - 4.0*Cs2_RZF*Can2_RZF)) )

      Cmax_ZF[i,j,k] = max(0, u_LZF+c_LZF, u_RZF+c_RZF)
      Cmin_ZF[i,j,k] = min(0, u_LZF-c_LZF, u_RZF-c_RZF)
   end

end


end
