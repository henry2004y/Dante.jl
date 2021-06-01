# Input Parameters

Here is a list of available input arguments in the TOML file.

## Parameters
* nD: dimension of the system
  * 1, 2, 3
* Scheme: numerical schemes to be used
  * "Rusanov"
  * "HLLE"
* Order: order of accuracy
  * 1
  * 2
* CFL: stability control parameter, (0,1)
* limiter: slope limiter for the 2nd order schemes
  * "MM": minmod
  * "MC":
* TimeAccurate: logical for running in time accurate mode
  * true
  * false
* UseConservative: logical for using energy conservation equation
  * true
  * false
* IC: initial conditions
  * "density wave"
  * "square wave"
  * "contact discontinuity"
  * "Riemann"
* RiemannProblemType: [1,12]
* nStep: total number of steps
* tEnd: end time in time accurate mode

## Grid
* TypeGrid: coordinate system
  * "Cartesian"
* xyzMinMax: range of the coordinates
  * [[0.0, 1.0]]
* nI: number of cells in the first dimension
* nJ: number of cells in the second dimension
* nK: number of cells in the third dimension
* BCtype: boundary conditions
  * ["float", "float"]
  * ["periodic", "periodic"]

## Plots
* DoPlot: logicals for plotting
  * true
  * false
* PlotVar: variable name to be plotted
  * "rho"
  * "ux"
  * "p"
* PlotInterval: plotting frequency in steps
* PlotType: type of plots
  * ["x", "1D"]
