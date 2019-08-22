module Divergence

export divergence_ndgrid

function divergence_ndgrid(hx::StepRangeLen,hy::StepRangeLen,hz::StepRangeLen,
   var::Array{Float64,4})::Array{Float64,3}

   siz = size(var)::NTuple{4,Int64}

   px = zeros(Float64,Base.front(siz))
   qy = zeros(Float64,Base.front(siz))
   rz = zeros(Float64,Base.front(siz))

   n = size(hx,1)
   # Right now do nothing for the ghost cells maybe needed later!
   # Take central differences on interior points
   if n > 2
      compute_px(siz,px,n,hx,var)
   end

   n = size(hy,1)
   if n > 2
      compute_qy(siz,qy,n,hx,var)
   end

   n = size(hz,1)
   if n > 2
      compute_rz(siz,rz,n,hx,var)
   end

   div = px .+ qy .+ rz
end

function compute_px(siz,px,n,hx,var)
   @inbounds for i=1:n-2, j=1:siz[2], k=1:siz[3]
         px[i+1,j,k] = (var[i+2,j,k,1] - var[i,j,k,1])/(hx[i+2] - hx[i])
   end
end

function compute_qy(siz,qy,n,hy,var)
   @inbounds for i=1:siz[1], j=1:n-2, k=1:siz[3]
         qy[i,j+1,k] = (var[i,j+1,k,2] - var[i,j,k,2])/(hy[i+2] - hy[i])
   end
end

function compute_rz(siz,rz,n,hz,var)
   @inbounds for i=1:siz[1], j=1:siz[2], k=1:n-2
         rz[i,j,k+1] = (var[i,j,k+1,3] - var[i,j,k,3])/(hz[i+2] - hz[i])
   end
end

end
