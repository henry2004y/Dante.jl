var documenterSearchIndex = {"docs":
[{"location":"input/#Input-Parameters","page":"Input Parameters","title":"Input Parameters","text":"","category":"section"},{"location":"input/","page":"Input Parameters","title":"Input Parameters","text":"Here is a list of available input arguments in the TOML file.","category":"page"},{"location":"input/#Parameters","page":"Input Parameters","title":"Parameters","text":"","category":"section"},{"location":"input/","page":"Input Parameters","title":"Input Parameters","text":"nD: dimension of the system\n1, 2, 3\nScheme: numerical schemes to be used\n\"Rusanov\"\n\"HLLE\"\nOrder: order of accuracy\n1\n2\nCFL: stability control parameter, (0,1)\nlimiter: slope limiter for the 2nd order schemes\n\"MM\": minmod\n\"MC\": monotonized central difference\nTimeAccurate: logical for running in time accurate mode\ntrue\nfalse\nUseConservative: logical for using energy conservation equation\ntrue\nfalse\nIC: initial conditions\n\"density wave\"\n\"square wave\"\n\"contact discontinuity\"\n\"Riemann\"\nRiemannProblemType: [1,12]\nnStep: total number of steps\ntEnd: end time in time accurate mode","category":"page"},{"location":"input/#Grid","page":"Input Parameters","title":"Grid","text":"","category":"section"},{"location":"input/","page":"Input Parameters","title":"Input Parameters","text":"TypeGrid: coordinate system\n\"Cartesian\"\nxyzMinMax: range of the coordinates\n[[0.0, 1.0]]\nnI: number of cells in the first dimension\nnJ: number of cells in the second dimension\nnK: number of cells in the third dimension\nBCtype: boundary conditions\n[\"float\", \"float\"]\n[\"periodic\", \"periodic\"]","category":"page"},{"location":"input/#Plots","page":"Input Parameters","title":"Plots","text":"","category":"section"},{"location":"input/","page":"Input Parameters","title":"Input Parameters","text":"DoPlot: logicals for plotting\ntrue\nfalse\nPlotVar: variable name to be plotted\n\"rho\"\n\"ux\"\n\"p\"\nPlotInterval: plotting frequency in steps\nPlotType: type of plots\n[\"x\", \"1D\"]","category":"page"},{"location":"api/#API","page":"API Reference","title":"API","text":"","category":"section"},{"location":"api/#Public","page":"API Reference","title":"Public","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [Dante]\nPrivate = false\nOrder = [:constant, :type, :function]","category":"page"},{"location":"api/#Dante.setInitRiemann","page":"API Reference","title":"Dante.setInitRiemann","text":"setInitRiemann(RiemannProblemType, Verbose)\n\nSet the initial conditions of Riemann problems. Note that currently tEnd and CFL number can only be set in PARAM.toml, and cannot be changed afterwards.\n\n\n\n\n\n","category":"function"},{"location":"api/#Dante.solve","page":"API Reference","title":"Dante.solve","text":"solve(paramFile)\n\nRun the model given input parameter file paramFile.\n\n\n\n\n\n","category":"function"},{"location":"api/#Private","page":"API Reference","title":"Private","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [Dante]\nPublic = false\nOrder = [:constant, :type, :function]","category":"page"},{"location":"api/#Dante.Param3D","page":"API Reference","title":"Dante.Param3D","text":"Model input parameters\n\n\n\n\n\n","category":"type"},{"location":"api/#Dante.advance!-Tuple{Any, Any}","page":"API Reference","title":"Dante.advance!","text":"Explicit time advance.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.calc_face_value!-Tuple{Dante.Param, Any, Dante.FaceState, Dante.FaceGradient}","page":"API Reference","title":"Dante.calc_face_value!","text":"calc_face_value!(param, state_GV, faceState, faceGradient)\n\nType instability for the return type is introduced, but it seems ok. I don't know how to modify faceState for views in 1st order.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.divergence!-Tuple{Any, Any, Any}","page":"API Reference","title":"Dante.divergence!","text":"divergence!(param, vec, div)\n\nCalculate the divergence of vectors specialized to my grid size. Always assume starting with i -> j -> k for 1/2/3D! Right now do nothing for the ghost cells. Maybe needed later! Take central differences on interior points.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.divergence_ndgrid!-NTuple{5, Any}","page":"API Reference","title":"Dante.divergence_ndgrid!","text":"divergence_ndgrid!(hx, hy, hz, vec, div)\n\nGeneric divergence div of vector vec with step lengths hx, hy, hz.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.getEulerExactSol-NTuple{8, Any}","page":"API Reference","title":"Dante.getEulerExactSol","text":"getEulerExactSol(ρ1, u1, p1, ρ4, u4, p4, tEnd, n)\n\nClassical Gas Exact Riemann Solver for shock tube problems. This programs is based on the code of Principles Of Computational Fluid Dynamics by P. Wesseling.\n\nnote: Note\nA Cavitation Check is incorporated in the code. It further prevents plotting for possible but physically unlikely case of expansion shocks.\n\nArguments\n\nρ1::Float64: left density at t=0.\nu1::Float64: left velocity at t=0.\np1::Float64: left pressure at t=0.\nρ4::Float64: right density at t=0.\nu4::Float64: right velocity at t=0.\np4::Float64: right pressure at t=0.\ntEnd::Float64: final solution time.\nn::Integer64: the gas Degree of freedom.\n\nCoded by Manuel Diaz, IAM, NTU 03/09/2011. Migrated by Hongyang Zhou from MATLAB to Julia, 11/05/2019\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.get_speed_max!-Tuple{Dante.Param, Dante.FaceState, Dante.SpeedFlux}","page":"API Reference","title":"Dante.get_speed_max!","text":"Calculate the maximum speed in each direction.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.minmod-Tuple{Any, Any, Any}","page":"API Reference","title":"Dante.minmod","text":"minmod(a, b, c)\n\nFor three inputs, use Harten's generalized definition. Return zero if opposite sign, otherwise the one of smaller magnitude.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.minmod-Tuple{Any, Any}","page":"API Reference","title":"Dante.minmod","text":"minmod(a, b)\n\nReturn zero if opposite sign, otherwise the one of smaller magnitude.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.plot_Riemann_exact-Tuple{Dante.Param}","page":"API Reference","title":"Dante.plot_Riemann_exact","text":"Plot the analytical solution of shock tube problem.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.plotvar-Tuple{Dante.Param, Any, Any}","page":"API Reference","title":"Dante.plotvar","text":"Plot 1D variables along a line.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.setParameters-Tuple{Any}","page":"API Reference","title":"Dante.setParameters","text":"setParameters(filename)\n\nRead parameters from PARAM.toml and construct the parameter list.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.update_state!-Tuple{Dante.Param, Any, AbstractFloat, Dante.FaceFlux, Any}","page":"API Reference","title":"Dante.update_state!","text":"Time-accurate state update.\n\n\n\n\n\n","category":"method"},{"location":"api/#Dante.update_state!-Tuple{Dante.Param, Any, Any, Dante.FaceFlux, Any}","page":"API Reference","title":"Dante.update_state!","text":"Local timestepping state update.\n\n\n\n\n\n","category":"method"},{"location":"getting_started/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"To install the package,","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"using Pkg; Pkg.add(PackageSpec(url=\"https://github.com/henry2004y/Dante\", rev=\"master\"))","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"To run with a given set of parameters,","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"solve(PARAMfile)","category":"page"},{"location":"fv/#Finite-Volume-Method","page":"Finite Volume","title":"Finite Volume Method","text":"","category":"section"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"The numerical approximation comes in when we define state variables at cell centers and fluxes at cell faces: the question of how to get the flux on faces from cell-centered states is the key to finite volumn method. In BATS-R-US, we use blocks to organize cells and distribute them among MPI processes.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"(Image: Figure 1)","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"shows the simplest block structure in 3D Cartesian coordinates consisting of 10times10times10 cells.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"(Image: Figure 2)","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"shows one cell in 3D Cartesian coordinates. For each face of the cell, we use left and right states to denote the states on the two sides of interface.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"Let us denote the approximate fluxes at faces with widehatF, the left states with subscript L and right states with subscript R. The simplest flux function in MHD is the local Lax-Friedrichs, or Rusanov flux, given by","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"widehatmathbfF = frac12(mathbfF_L+mathbfF_R)-frac12 c_max(mathbfU_R-mathbfU_L)","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"where mathbfF_L = mathbfF(mathbfU_L) mathbfF_R = mathbfF(mathbfU_R), and c_max is the maximum speed of any wave in the system, considering both mathbfU_L and mathbfU_R. (It is the same for all flux calculations on this face, but different at different faces!) In ideal MHD, this maximum wave speed can be expressed as","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"c_max =c_max(mathbfU_LmathbfU_Rwidehatn_ie) =  bar u_n + sqrt fracc_f^2+  sqrtc_f^4 - 4 c_s^2 c^2_An2","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"where the fast magnetosonic speed squared is","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"  c_f^2 = c_s^2 + c_A^2","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"the sound speed squared is","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"  c_s^2 = fracgamma bar pbar rho=fracgamma (p_L+p_R)rho_L+rho_R","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"the Alfvén speed squared is","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"  c_A^2 = fracbar B^2bar rho = frac(B_xL+B_xR)^2+(B_yL+B_yR)^2+(B_zL+B_zR)^22(rho_L+rho_R)","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"and the normal component of the Alfvén speed squared is ","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"  c_An^2 = fracbar B_n^2bar rho = c_A^2 widehatn_f","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"note: Note\nPhysically the Alfvén velocity is written as mathbfc_A = fracmathbfB^2rho(mu_0 is absorbed into mathbfB in the normalized units). At cell faces, we take the average of the magnetic field and densities between neighboring cell centers to get the velocity on the faces: mathbfc_textAface = fracbarmathbfB^2barrho The normal component of the Alfvén velocity to the face is then mathbfc_An =  frac(barmathbfBcdotwidehatn_f)^2barrho and Alfvén speed is just the magnitude of this.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"In the discretization form, assuming we are using first-order forward Euler method in time, integrating over the whole cell i and taking cell averages at the center, the compact conservative form of the MHD equations gives","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"V_ifrac1Delta t_i^n(u_i^n+1-u_i^n)+sum_f=1^6widehatF(u_iu_ewidehatn_f)A_f = V_i S_i^n","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"where V_i is the volume of the i^th cell, A_f is the face area of f (between cell i and its adjacent cell e), widehatF is the flux function between cell i and cell e, u_f is the state variables in adjacent cell e, widehatn_f is the face normal vector pointing from cell i to cell e, and superscript n is the timestep. The second term on the LHS gives the integrated normal flux of state u on face between cell i and cell e. Rewrite the above equation, we get","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"u_i^n+1  = u_i^n - fracDelta t_i^nV_iBig( sum_f=1^6widehatF(u_iu_ewidehatn_f)A_fBig) + Delta t_i^n S_i^n ","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"math which is the explicit update equation for state variables from timestep n to n+1.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"Finally, to set the discrete timestep, we apply the CFL condition for numerical stability. ","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"Delta t_i = textCFLfracV_is_ixA_ix + s_iyA_iy + s_izA_iz= textCFL frac1fracs_ixDelta x+fracs_iyDelta y+fracs_izDelta z","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"where CFL number is a constant over all cells, prescribed by the user, and s_i is the maximum wave speed in cell i. For general cases, this is the same as c_max used for flux function.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"warning: Warning\nBe careful: this time the maximum wave speed is defined in cell centers, not faces. So instead of using the average states between two cells, you only need cell center states of i to compute this s_i. In my first MATLAB version of the code, I made a mistake initially by using the maximum of the face c_max instead of calculating cell-centered c_max directly. This resulted in oscillation behaviors in certain regions.","category":"page"},{"location":"fv/","page":"Finite Volume","title":"Finite Volume","text":"note: Note\nDante is a simplified version of BATS-R-US, which is a generalized code that can handle different coordinates with different mesh sizes. Therefore we shouldn't simplify the calculations too much even though it is good as a start. Also note that what`s been shown here is the simplest first-order explicit scheme: we have much more complicated schemes with higher-order and a mixture of implicit-explicit implementations.","category":"page"},{"location":"test/#Standard-Test","page":"Standard Test","title":"Standard Test","text":"","category":"section"},{"location":"test/","page":"Standard Test","title":"Standard Test","text":"Brio-Wu shocktube test is an extension of Sods shock test in hydrodynamics. It adds magnetic fields to the fluid system.","category":"page"},{"location":"#Dante.jl","page":"Home","title":"Dante.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Dante.jl is a toy finite volume ideal MHD model on structured mesh. It is the descendant of the original version in MATLAB.","category":"page"},{"location":"#Benchmark","page":"Home","title":"Benchmark","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Riemann test 1 (Sods problem) without plotting (test/PARAM_test.toml):","category":"page"},{"location":"","page":"Home","title":"Home","text":"4.859 ms (1182 allocations: 691.80 KiB)","category":"page"},{"location":"mhd/#Magnetohydrodynamics-Model","page":"MHD","title":"Magnetohydrodynamics Model","text":"","category":"section"},{"location":"mhd/","page":"MHD","title":"MHD","text":"note: Note\nAs of June 2021, KaTeX lacks the full support for equation numbering.","category":"page"},{"location":"mhd/#Ideal-MHD","page":"MHD","title":"Ideal MHD","text":"","category":"section"},{"location":"mhd/","page":"MHD","title":"MHD","text":"The ideal MHD equations can be written in (near) conservative form as","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"beginaligned\nfracpartial rhopartial t + nablacdot(rho mathbfu) = 0 \nfracpartial rho mathbfupartial t + nablacdot left rho mathbfu mathbfu + barbarI (p + frac12 B^2)\n- mathbfBmathbfB right = -mathbfBnablacdotmathbfB \nfracpartial mathbfBpartial t + nablacdot (mathbfumathbfB - mathbfBmathbfu)\n= -mathbfunablacdotmathbfB \nfracpartial epartial t + nablacdotleft mathbfu ( e + p + frac12 B^2 )\n-mathbfucdotmathbfBmathbfB right = -mathbfucdotmathbfBnablacdotmathbfB\nendaligned","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"where rho is mass density, mathbfu is velocity, p is pressure, mathbfB is the magnetic field and barbarI is the identity matrix.  The total energy density is","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"  tagtotalE\n  e = fracpgamma-1 + fracrho u^22 + fracB^22","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"where gamma is the adiabatic index.","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"note: Note\nWe keep the divergence of the magnetic field on the RHS because of numerical accuracy. Even though physically valid by far in nature, it is not guaranteed to be zero in a numerical model.","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"We can also substitute the energy equation with pressure equation for a non-conservative form:","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"  tagpressure\n  fracpartial ppartial t + nablacdot (p mathbfu) = -(gamma - 1) p nabla cdot mathbfu ","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"The source terms proportional to nablacdotmathbfB and nablacdotmathbfu on the right hand sides are not evaluated as a flux, but the ingredients (normal components of mathbfB and mathbfu are calculated and used later. In total, there are 8 primitive variables (rhomathbfumathbfbp). Note that in BAT-S-RUS, the momentum rhomathbfu is actually stored instead of velocity mathbfu. ","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"In Cartesian coordinates, we can express the above equations in more detailed and compact form. If we use mathbfA to represent all the state variables, we have","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"fracpartialmathbfApartial t+nablacdotvecmathbfF(mathbfB)=mathbfS","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"by defining the state and flux vectors as","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"beginaligned\ntextstate  mathbfA = \n\tbeginbmatrix\n\trho  rho u_x  rho u_y  rho u_z  B_x  B_y  B_z  p  e \n\tendbmatrix\nquad\ntextflux  vecmathbfF = \n\tbeginbmatrix\n\trho u_x  rho u_x u_x + (p+frac12B^2)- B_xB_x  rho u_x u_y - B_x B_y  rho u_x u_z - B_xB_z  0  u_x B_y - B_x u_y \n\tu_x B_z - B_x u_z  p u_x u_x (e+p+frac12B^2) - (u_xB_x+u_yB_y+u_zB_z)B_x\n\tendbmatrixwidehatx + \n\tbeginbmatrix\n\trho u_y  rho u_y u_x - B_yB_x  rho u_y u_y + (p+frac12B^2)- B_y B_y  rho u_y u_z - B_yB_z  u_y B_x - B_y u_x  0 \n\tu_y B_z - B_y u_z  p u_y  u_y (e+p+frac12B^2) - (u_xB_x+u_yB_y+u_zB_z)B_y\n\tendbmatrixwidehaty +\n\tbeginbmatrix\n\trho u_z  rho u_z u_x - B_zB_x  rho u_z u_y - B_z B_y  rho u_z u_z + (p+frac12B^2) - B_zB_z  u_z B_x - B_z u_x   u_z B_y - B_z u_y  0  p u_z u_z (e+p+frac12B^2) - (u_xB_x+u_yB_y+u_zB_z)B_z\n\tendbmatrixwidehatz \ntextsource  mathbfS = \n\tbeginbmatrix\n\t0  \n\t-B_x(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\t-B_y(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\t-B_z(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\t-u_x(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\t-u_y(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\t-u_z(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\t-(gamma-1)p(fracpartial u_xpartial x+fracpartial u_ypartial y+fracpartial u_zpartial z) \n\t-(u_xB_x+u_yB_y+u_zB_z)(fracpartial B_xpartial x+fracpartial B_ypartial y+fracpartial B_zpartial z) \n\tendbmatrix\nendaligned","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"Note that the last two equations overspecify the system. We only need to choose one of them.","category":"page"},{"location":"mhd/","page":"MHD","title":"MHD","text":"Up to this stage, all the equations above are exact. The next question is how do we solve this system numerically?","category":"page"}]
}