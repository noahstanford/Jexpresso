using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
import SciMLBase

include("../abstractTypes.jl")
include("../../io/write_output.jl")

function time_loop!(SD,
                    QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M,
                    De, Le,
                    Nt, Δt,
                    neqns, 
                    inputs::Dict,
                    BCT,
                    OUTPUT_DIR::String,
                    T)
    
    #
    # ODE
    #
    u      = zeros(T, mesh.npoin);
    u     .= qp.qn[:,1];
    tspan  = (inputs[:tinit], inputs[:tend])    
    params = (; T, SD, QT, PT, BCT, neqns, basis, ω, mesh, metrics, inputs, M, De, Le)
    prob   = ODEProblem(rhs!,
                        u,
                        tspan,
                        params);
    
    
    println(" # Solving ODE with ................................\n")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    
    @time    sol = solve(prob,
                         inputs[:ode_solver],
                         dt = Δt,
                         saveat = range(T(0.), Nt*T(Δt), length=inputs[:ndiagnostics_outputs]),
                         progress = true,
                         progress_message = (dt, u, p, t) -> t)
    println(" # Solving ODE with  ................................ DONE\n")
    
    write_output(sol, SD, mesh, OUTPUT_DIR, inputs; outformat=inputs[:outformat])
    
    
end
