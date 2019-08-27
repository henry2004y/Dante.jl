using BenchmarkTools

function compute_px1(siz,px,n,hx,var)
   for k=1:siz[3], j=1:siz[2], i=1:n-2
         px[i+1,j,k] = (var[i+2,j,k,1] - var[i,j,k,1])/(hx[i+2] - hx[i])
   end
end
function compute_px2(siz,px,n,hx,var)
   for i=1:n-2, j=1:siz[2], k=1:siz[3]
         px[i+1,j,k] = (var[i+2,j,k,1] - var[i,j,k,1])/(hx[i+2] - hx[i])
   end
end

function compute_px3(siz,px,n,hx,var)
   @inbounds for k=1:siz[3], j=1:siz[2], i=1:n-2
         px[i+1,j,k] = (var[i+2,j,k,1] - var[i,j,k,1])/(hx[i+2] - hx[i])
   end
end
function compute_px4(siz,px,n,hx,var)
   @inbounds for i=1:n-2, j=1:siz[2], k=1:siz[3]
         px[i+1,j,k] = (var[i+2,j,k,1] - var[i,j,k,1])/(hx[i+2] - hx[i])
   end
end
compute_px5(siz,px,n,hx,var) = @. px[2:n-1,:,:] = (var[3:n,:,:,1] - var[1:n-2,:,:,1])/(hx[3:n] - hx[1:n-2])
compute_px6(siz,px,n,hx,var) = @inbounds @. px[2:n-1,:,:] = (var[3:n,:,:,1] - var[1:n-2,:,:,1])/(hx[3:n] - hx[1:n-2])
compute_px7(siz,px,n,hx,var) = @views @. px[2:n-1,:,:] = (var[3:n,:,:,1] - var[1:n-2,:,:,1])/(hx[3:n] - hx[1:n-2])


function foo(hx::StepRangeLen,hy::StepRangeLen,hz::StepRangeLen,
   var::Array{Float64,4})::Array{Float64,3}

   siz = size(var)::NTuple{4,Int64}

   px = zeros(Float64,siz[1:3]...)
   qy = zeros(Float64,siz[1:3]...)
   rz = zeros(Float64,siz[1:3]...)

   n = size(hx,1)

   if n > 2
      @btime compute_px1($siz,$px,$n,$hx,$var)
      @btime compute_px2($siz,$px,$n,$hx,$var)
      @btime compute_px3($siz,$px,$n,$hx,$var)
      @btime compute_px4($siz,$px,$n,$hx,$var)
      @btime compute_px5($siz,$px,$n,$hx,$var)
      @btime compute_px6($siz,$px,$n,$hx,$var)
      @btime compute_px7($siz,$px,$n,$hx,$var)
      @time compute_px1(siz,px,n,hx,var)
      @time compute_px2(siz,px,n,hx,var)
      @time compute_px3(siz,px,n,hx,var)
      @time compute_px4(siz,px,n,hx,var)
      @time compute_px5(siz,px,n,hx,var)
      @time compute_px6(siz,px,n,hx,var)
      @time compute_px7(siz,px,n,hx,var)
   end

   px
end

function foo_stable(hx::StepRangeLen,hy::StepRangeLen,hz::StepRangeLen,
   var::Array{Float64,4})::Array{Float64,3}

   siz = size(var)::NTuple{4,Int64}

   px = zeros(Float64,Base.front(siz))
   qy = zeros(Float64,Base.front(siz))
   rz = zeros(Float64,Base.front(siz))

   n = size(hx,1)

   if n > 2
      @btime compute_px1($siz,$px,$n,$hx,$var)
      @btime compute_px2($siz,$px,$n,$hx,$var)
      @btime compute_px3($siz,$px,$n,$hx,$var)
      @btime compute_px4($siz,$px,$n,$hx,$var)
      @btime compute_px5($siz,$px,$n,$hx,$var)
      @btime compute_px6($siz,$px,$n,$hx,$var)
      @btime compute_px7($siz,$px,$n,$hx,$var)
      @time compute_px1(siz,px,n,hx,var)
      @time compute_px2(siz,px,n,hx,var)
      @time compute_px3(siz,px,n,hx,var)
      @time compute_px4(siz,px,n,hx,var)
      @time compute_px5(siz,px,n,hx,var)
      @time compute_px6(siz,px,n,hx,var)
      @time compute_px7(siz,px,n,hx,var)
   end

   px
end

function test()
   x = range(0.0, 1.0, length=10000)
   y = 0.5:0.5:1.5
   z = 0.5:0.5:1.5
   var = zeros(10000,3,3,3)
   px = foo(x,y,z,var);
   px = foo_stable(x,y,z,var);
end
test()
