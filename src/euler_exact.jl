using Roots

"""
	EulerExact(ρ1, u1, p1, ρ4, u4, p4, tEnd, n)

Classical Gas Exact Riemann Solver for shock tube problems
This programs is based on the code of
P. Wesseling. PRINCIPLES OF COMPUTATIONAL FLUID DYNAMICS Springer-Verlag, Berlin etc., 2001.
ISBN 3-540-67853-0
NOTE:
A Cavitation Check is incorporated in the code. It further prevents plotting for
 possible but physically unlikely case of expansion shocks.
# Arguments
- `ρ1::Float64`: left density at t=0.
- `u1::Float64`: left velocity at t=0.
- `p1::Float64`: left pressure at t=0.
- `ρ4::Float64`: right density at t=0.
- `u4::Float64`: right velocity at t=0.
- `p4::Float64`: right pressure at t=0.
- `tEnd::Float64`: final solution time.
- `n::Integer64`: the gas Degree of freedom.
Coded by Manuel Diaz, IAM, NTU 03/09/2011.
Migrated by Hongyang Zhou from MATLAB to Julia, 11/05/2019
"""
function EulerExact(ρ1, u1, p1, ρ4, u4, p4, tEnd, n)

   # Gamma values
   γ = (n+2)/n; α = (γ+1)/(γ-1)

   # Assumed structure of exact solution
   #
   #    \         /      |con |       |s|
   #     \   f   /       |tact|       |h|
   # left \  a  /  state |disc| state |o| right
   # state \ n /    2    |cont|   3   |c| state
   #   1    \ /          |tinu|       |k|   4
   #         |           |ity |       | |

   PRL      = p4/p1
   cright   = √(γ*p4/ρ4)
   cleft    = √(γ*p1/ρ1)
   CRL      = cright/cleft
   MachLeft = (u1-u4)/cleft

   # Basic shock tube relation equation (10.51)
   f(P) = (1 + MachLeft*(γ-1)/2 - (γ-1)*CRL*(P-1)/√(2*γ*(γ-1+(γ+1)*P)))^(2*γ/
      (γ-1))/P - PRL

   # solve for P = p34 = p3/p4
   p34 = find_zero(f, 3)

   p3 = p34*p4
   ρ3 = ρ4*(1+α*p34)/(α+p34)
   ρ2 = ρ1*(p34*p4/p1)^(1/γ)
   u2 = u1 -u4 + (2/(γ-1))*cleft*(1-(p34*p4/p1)^((γ-1)/(2*γ)))
   c2 = √(γ*p3/ρ2)
   spos = 0.5 + tEnd*cright*√((γ-1)/(2*γ)+(γ+1)/(2*γ)*p34) + tEnd*u4

   x0 = 0.5
   conpos=x0 + u2*tEnd + tEnd*u4	   # Position of contact discontinuity
   pos1 = x0 + (u1 - cleft)*tEnd	   # Start of expansion fan
   pos2 = x0 + (u2 + u4 - c2)*tEnd	# End of expansion fan

   # Plot structures
   x = 0:0.002:1
   xsize = size(x)
   p = zeros(xsize)
   ux= zeros(xsize)
   ρ = zeros(xsize)
   Mach = zeros(xsize)
   cexact = zeros(xsize)

   for i = 1:length(x)
      if x[i] ≤ pos1
         ρ[i], p[i], ux[i] = ρ1, p1, u1
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      elseif x[i] ≤ pos2
         p[i] = p1*(1+(pos1-x[i])/(cleft*α*tEnd))^(2*γ/(γ-1))
         ρ[i] = ρ1*(1+(pos1-x[i])/(cleft*α*tEnd))^(2/(γ-1))
         ux[i] = u1 + (2/(γ+1))*(x[i]-pos1)/tEnd
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      elseif x[i] ≤ conpos
         ρ[i], p[i], ux[i] = ρ2, p3, u2+u4
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      elseif x[i] ≤ spos
         ρ[i], p[i], ux[i] = ρ3, p3, u2+u4
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      else
         ρ[i], p[i], ux[i] = ρ4, p4, u4
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      end
   end
   entro = @. log(p/ρ^γ) # entropy
   e = @. p/((γ-1)*ρ)	 # internal energy
   t = @. 2/n*e          # temperature

   return x, ρ, ux, p, e, t, Mach, entro
end

# Vector x version
function EulerExact(ρ1, u1, p1, ρ4, u4, p4, tEnd, n, x)

   # Gamma values
   γ=(n+2)/n; α=(γ+1)/(γ-1)

   # Assumed structure of exact solution
   #
   #    \         /      |con |       |s|
   #     \   f   /       |tact|       |h|
   # left \  a  /  state |disc| state |o| right
   # state \ n /    2    |cont|   3   |c| state
   #   1    \ /          |tinu|       |k|   4
   #         |           |ity |       | |

   PRL      = p4/p1
   cright   = √(γ*p4/ρ4)
   cleft    = √(γ*p1/ρ1)
   CRL      = cright/cleft
   MachLeft = (u1-u4)/cleft

   # Basic shock tube relation equation (10.51)
   f(P) = (1 + MachLeft*(γ-1)/2 - (γ-1)*CRL*(P-1)/√(2*γ*(γ-1+(γ+1)*P)))^(2*γ/
      (γ-1))/P - PRL

   # solve for P = p34 = p3/p4
   p34 = find_zero(f, 3)

   p3 = p34*p4
   ρ3 = ρ4*(1+α*p34)/(α+p34)
   ρ2 = ρ1*(p34*p4/p1)^(1/γ)
   u2 = u1 -u4 + (2/(γ-1))*cleft*(1-(p34*p4/p1)^((γ-1)/(2*γ)))
   c2 = √(γ*p3/ρ2)
   spos = 0.5 + tEnd*cright*√((γ-1)/(2*γ)+(γ+1)/(2*γ)*p34) + tEnd*u4

   x0 = 0.5
   conpos=x0 + u2*tEnd + tEnd*u4	   # Position of contact discontinuity
   pos1 = x0 + (u1 - cleft)*tEnd	   # Start of expansion fan
   pos2 = x0 + (u2 + u4 - c2)*tEnd	# End of expansion fan

   # Plot structures
   xsize = size(x)
   p = zeros(xsize)
   ux= zeros(xsize)
   ρ = zeros(xsize)
   Mach = zeros(xsize)
   cexact = zeros(xsize)

   for i = 1:length(x)
      if x[i] ≤ pos1
         ρ[i], p[i], ux[i] = ρ1, p1, u1
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      elseif x[i] ≤ pos2
         p[i] = p1*(1+(pos1-x[i])/(cleft*α*tEnd))^(2*γ/(γ-1))
         ρ[i] = ρ1*(1+(pos1-x[i])/(cleft*α*tEnd))^(2/(γ-1))
         ux[i] = u1 + (2/(γ+1))*(x[i]-pos1)/tEnd
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      elseif x[i] ≤ conpos
         ρ[i], p[i], ux[i] = ρ2, p3, u2+u4
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      elseif x[i] ≤ spos
         ρ[i], p[i], ux[i] = ρ3, p3, u2+u4
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      else
         ρ[i], p[i], ux[i] = ρ4, p4, u4
         cexact[i] = √(γ*p[i]/ρ[i])
         Mach[i] = ux[i]/cexact[i]
      end
   end
   entro = @. log(p/ρ^γ) # entropy
   e = @. p/((γ-1)*ρ)	 # internal energy
   t = @. 2/n*e          # temperature

   return x, ρ, ux, p, e, t, Mach, entro
end
