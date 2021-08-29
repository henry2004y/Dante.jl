# Input parameters

using TOML

abstract type Param end

const Rho_       = 1
const Ux_        = 2
const Uy_        = 3
const Uz_        = 4
const Mx_        = 2
const My_        = 3
const Mz_        = 4
const Bx_        = 5
const By_        = 6
const Bz_        = 7
const P_         = 8
const E_         = P_
const U_         = Ux_:Uz_
const M_         = Mx_:Mz_
const B_         = Bx_:Bz_

const γ = 5.0/3.0

"Model input parameters"
@with_kw struct Param3D{IntT, RealT, RangeT<:AbstractArray} <: Param
	nD::IntT
	nI::IntT
	nJ::IntT
	nK::IntT
	nG::IntT
	TypeGrid::String
	Order::Int = 2
	Scheme::String
	limiter::String = "MM"
	TimeAccurate::Bool = true
	UseConservative::Bool = true
	UseGPU::Bool = false
	verbose::Bool = false
	iMin::IntT
	iMax::IntT
	jMin::IntT
	jMax::IntT
	kMin::IntT
	kMax::IntT
	iMinAll::IntT
	iMaxAll::IntT
	jMinAll::IntT
	jMaxAll::IntT
	kMinAll::IntT
	kMaxAll::IntT
	CFL::RealT
	nStep::IntT             # Total timesteps
	tEnd::RealT
	TypeBc::Vector{String} # Type of boundary conditions
	IC::String
	RiemannProblemType::IntT

	nVar::IntT

	DoPlot::Bool = false
	PlotInterval::IntT = -1
	PlotVar::String
	x::RangeT
	y::RangeT
	z::RangeT
	FullSize::Vector{IntT}
	GridSize::Vector{IntT}
	CellSize_D::Vector{RealT}
end

"""
	setParameters(filename)

Read parameters from PARAM.toml and construct the parameter list.
"""
function setParameters(filename)

	paramIn = TOML.parsefile(filename)

	nD = paramIn["Parameters"]["nD"]::Int
	TypeGrid = paramIn["Grid"]["TypeGrid"]::String
	Order = paramIn["Parameters"]["Order"]::Int
	Scheme = paramIn["Parameters"]["Scheme"]::String
	limiter = paramIn["Parameters"]["limiter"]::String
	nI = paramIn["Grid"]["nI"]::Int
	nJ = paramIn["Grid"]["nJ"]::Int
	nK = paramIn["Grid"]["nK"]::Int
	xyzMinMax_D = hcat(paramIn["Grid"]["xyzMinMax"]...)::Array{Float64}

	CoordMinMax_D = [0.0 0.0 0.0; 1.0 1.0 1.0]

	for i in eachindex(xyzMinMax_D)
		CoordMinMax_D[i] = xyzMinMax_D[i]
	end

	CFL = paramIn["Parameters"]["CFL"]::Float64
	TimeAccurate = paramIn["Parameters"]["TimeAccurate"]::Bool
	UseConservative = paramIn["Parameters"]["UseConservative"]::Bool
	UseGPU = paramIn["Parameters"]["UseGPU"]::Bool
	nStep  = paramIn["Parameters"]["nStep"]::Int64
	tEnd = paramIn["Parameters"]["tEnd"]::Float64
	@assert(0 ≤ tEnd, "Simulation time must be positive!")
	DoPlot = paramIn["Plots"]["DoPlot"]::Bool
	PlotInterval = paramIn["Plots"]["PlotInterval"]::Int64
	PlotVar = paramIn["Plots"]["PlotVar"]::String
	IC = paramIn["Parameters"]["IC"]::String
	RiemannProblemType = paramIn["Parameters"]["RiemannProblemType"]::Int

	nG         = Order
	nStage     = Order

	@assert(1 ≤ nD ≤ 3, "Dimension error, nD = $(nD)!")
	if nD == 3
		GridSize = [nI, nJ, nK]
	elseif nD == 2
		GridSize = [nI, nJ, 1]
	elseif nD == 1
		GridSize = [nI, 1, 1]
	end
	if TypeGrid == "Cartesian"
		#X,Y,Z = ndgrid(x,y,z)
	elseif TypeGrid == "Spherical"

	else
		error("Unknown TypeGrid = $(TypeGrid)")
	end

	CellSize_D = (CoordMinMax_D[2,:] .- CoordMinMax_D[1,:]) ./ GridSize
	# Size including ghost cells
	FullSize   = GridSize .+ 2*nG

	x = range(CoordMinMax_D[1,1] + (0.5-nG)*CellSize_D[1],
	CoordMinMax_D[2,1] + (nG-0.5)*CellSize_D[1], length=FullSize[1])
	y = range(CoordMinMax_D[1,2] + (0.5-nG)*CellSize_D[2],
	CoordMinMax_D[2,2] + (nG-0.5)*CellSize_D[2], length=FullSize[2])
	z = range(CoordMinMax_D[1,3] + (0.5-nG)*CellSize_D[3],
	CoordMinMax_D[2,3] + (nG-0.5)*CellSize_D[3], length=FullSize[3])

	# Indexes for physical cells
	iMin       = 1 + nG
	iMax       = nI + nG
	# Indexes for all cells including ghost cells
	iMinAll    = 1
	iMaxAll    = nI + 2*nG

	jMin       = 1 + nG
	jMax       = nJ + nG
	jMinAll    = 1
	jMaxAll    = nJ + 2*nG

	kMin       = 1 + nG
	kMax       = nK + nG
	kMinAll    = 1
	kMaxAll    = nK + 2*nG

	nVar = 8
	BCtype = paramIn["Grid"]["BCtype"]::Vector{String}
	# This is needed for correctly updating ghost cell values in Boundary.
	if length(BCtype) < 6
		for i = 1:6-length(BCtype)
			push!(BCtype,"float")
		end
	end

	verbose = paramIn["Parameters"]["verbose"]::Bool

	param = Param3D(nD,nI,nJ,nK,nG,TypeGrid,Order,Scheme,limiter,TimeAccurate,
		UseConservative,UseGPU,verbose,
		iMin,iMax,jMin,jMax,kMin,kMax,iMinAll,iMaxAll,jMinAll, jMaxAll,kMinAll,kMaxAll,
		CFL,nStep,tEnd,BCtype,IC,RiemannProblemType,nVar,DoPlot,
		PlotInterval,PlotVar,x,y,z,FullSize,GridSize,CellSize_D)
end