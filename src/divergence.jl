# Divergence calculation

"""
	divergence_ndgrid!(hx, hy, hz, vec, div)

Generic divergence `div` of vector `vec` with step lengths `hx`, `hy`, `hz`.
"""
function divergence_ndgrid!(hx, hy, hz, vec, div)

   siz = size(vec)::NTuple{4,Int64}

   px = zeros(Float64, Base.front(siz))
   qy = zeros(Float64, Base.front(siz))
   rz = zeros(Float64, Base.front(siz))

   n = size(hx, 1)
   # Right now do nothing for the ghost cells maybe needed later!
   # Take central differences on interior points
   if n > 2
      compute_px(siz, px, n, hx, vec)
   end

   n = size(hy, 1)
   if n > 2
      compute_qy(siz, qy, n, hx, vec)
   end

   n = size(hz, 1)
   if n > 2
      compute_rz(siz, rz, n, hx, vec)
   end

   @. div = px + qy + rz
   return
end

function divergence_ndgrid(hx, hy, hz, vec)

   siz = size(vec)

   div = zeros(Float64, Base.front(siz))

   divergence_ndgrid!(hx, hy, hz, vec, div)

   return div
end

function compute_px(siz, px, n, hx, vec)
   @inbounds for i = 1:n-2, j = 1:siz[2], k = 1:siz[3]
      px[i+1,j,k] = (vec[i+2,j,k,1] - vec[i,j,k,1])/(hx[i+2] - hx[i])
   end
end

function compute_qy(siz, qy, n, hy, vec)
   @inbounds for i = 1:siz[1], j = 1:n-2, k = 1:siz[3]
      qy[i,j+1,k] = (vec[i,j+2,k,2] - vec[i,j,k,2])/(hy[i+2] - hy[i])
   end
end

function compute_rz(siz, rz, n, hz, vec)
   @inbounds for i = 1:siz[1], j = 1:siz[2], k = 1:n-2
      rz[i,j,k+1] = (vec[i,j,k+2,3] - vec[i,j,k,3])/(hz[i+2] - hz[i])
   end
end

"""
   divergence!(param, vec, div)

Calculate the divergence of vectors specialized to my grid size.
Always assume starting with i -> j -> k for 1/2/3D!
Right now do nothing for the ghost cells. Maybe needed later!
Take central differences on interior points.
"""
@inline function divergence!(param, vec, div)

   x, y, z = param.x, param.y, param.z
   nI, nJ, nK, nG = param.nI, param.nJ, param.nK, param.nG
   iMin, iMax, jMin, jMax, kMin, kMax =
   param.iMin, param.iMax, param.jMin, param.jMax, param.kMin, param.kMax
   
   div .= 0.0

   @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
      div[i-nG,j-nG,k-nG] = (vec[i+1,j,k,1] - vec[i-1,j,k,1]) / (x[i+1] - x[i-1])
   end

   if nJ > 1
      @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
         div[i-nG,j-nG,k-nG] += (vec[i,j+1,k,2] - vec[i,j-1,k,2]) / (y[i+1] - y[i-1])
      end
   end

   if nK > 1
      @inbounds for k = kMin:kMax, j = jMin:jMax, i = iMin:iMax
         div[i-nG,j-nG,k-nG] += (vec[i,j,k+1,3] - vec[i,j,k-1,3]) / (z[i+1] - z[i-1])
      end
   end

   return
end