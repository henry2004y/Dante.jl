# Test of Dante
using Dante, LinearAlgebra, Test

@testset "Riemann Problem" begin
   param, state_GV = main("PARAM_test.toml")
   # Obtain the initial states
   Rho, U, P, tEnd = set_init_Riemann(param.RiemannProblemType, false)
   # Exact solution
   xe, re, ue, pe, ee, te, Me, se =
      EulerExact(Rho[1], U[1], P[1], Rho[end], U[end], P[end], tEnd, 3,
      param.x[3:end-2])
   ϵ = norm(state_GV[3:end-2,3,3,8] .- pe)
   @show ϵ
   @test ϵ < 0.25
end
