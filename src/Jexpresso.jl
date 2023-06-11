module Jexpresso

using Revise
using BenchmarkTools
using Dates
using DelimitedFiles
using DataStructures
using UnicodePlots
using Printf
using Revise
using ElasticArrays
using StaticArrays
using StaticArrays: SVector, MVector, MArray, SMatrix, @SMatrix
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
using SnoopCompile
using SciMLBase: CallbackSet, DiscreteCallback,
                 ODEProblem, ODESolution, ODEFunction,
                 SplitODEProblem
import SciMLBase: get_du, get_tmp_cache, u_modified!,
                  AbstractODEIntegrator, init, step!, check_error,
                  get_proposed_dt, set_proposed_dt!,
                  terminate!, remake



const TInt   = Int64
const TFloat = Float64

#using DocStringExtensions

include(joinpath("equations", "AbstractEquations.jl"))

include(joinpath("kernel", "abstractTypes.jl"))

include(joinpath("kernel", "globalStructs.jl"))

include(joinpath("kernel", "physics", "globalConstantsPhysics.jl"))

include(joinpath("kernel", "physics", "constitutiveLaw.jl"))

include(joinpath("kernel", "infrastructure", "sem_setup.jl"))

include(joinpath("kernel", "boundaryconditions", "BCs.jl"))

include(joinpath("kernel", "operators", "rhs.jl"))

include(joinpath("kernel", "operators", "rhs_diff.jl"))

include(joinpath("kernel", "solvers", "TimeIntegrators.jl"))

include(joinpath("kernel", "solvers", "Axb.jl"))

include(joinpath("io", "mod_inputs.jl"))

include(joinpath("io", "write_output.jl"))

include("./run.jl")

end
