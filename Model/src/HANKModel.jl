
module HANKModel

export SSminstruct, SSstruct, Polstruct, PolEstruct, MITstruct,
       CalibTargets, SteadyStats, EntreTemp, Temp

export main_fct, options, parameters, grids, fractional_grid_index, steady_init, get_dt, find_z_stat, find_y_stat, get_y_states, exostate_transit_mat, exostate_transit_mat_entre, lab_fun, lev_fun, linear_solve!, HHsolve_entre!, HHsolve!, init_A, transit_mat_entre!, transit_mat!, iterHJBentre!, iterHJB!, fdepodrift!, depodrift!, upwind!, Aentreassemble!, Aassemble!, Bassemble!, minusATassemble!, bbuild!, sfbuild!, scebuild!, sdebuild!, cebuild!, findsteady, lab0update!, MIT_init, solveMIT!, solveMIT, MITupdate!, MITall!, initial_q_change!, change_distr, calibrate, calibobjective, calibstats, plotMIT!, plotMITconv, shock, simulate_states, debackersim, debackersim2, makedebacker, debacker_mobility_stats, percshare, stats, statsMIT, partial, bequestregimes!, findsteady_fixedB, beqStats, @set!, calibstats!, worker_nonresponse_by_joint_quantile, solve_phillips_curve_grid

export @set!, @load, @save, @unpack, load

using DelimitedFiles, SparseArrays, Distributions, LinearAlgebra, Random, Roots, BenchmarkTools, TimerOutputs, ILUZero, IterativeSolvers, Krylov, Setfield, Optim, MAT, Parameters, JuMP, Interpolations, JLD2, MultistartOptimization, Distributed, SharedArrays, CSV, Tables, BasicInterpolators, DataFrames, NLopt, Evolutionary, PrettyTables, BlackBoxOptim, NOMAD, ConjugateGradients
using Plots;
gr();

include("./MyGlobalOptim/MyGlobalOptim.jl")
using .MyGlobalOptim

include("ModelTypes.jl")
include("linear_solve!.jl")
include("Methods.jl")
using .Methods
include("options.jl")
include("parameters.jl")
include("sparse_utils.jl")
include("grids.jl")
include("steady_init.jl")
include("get_dt.jl")
include("find_z_stat.jl")
include("find_y_stat.jl")
include("get_y_states.jl")
include("exostate_transit_mat.jl")
include("exostate_transit_mat_entre.jl")
include("lab_fun.jl")
include("lev_fun.jl")
include("HHsolve_entre!.jl")
include("HHsolve!.jl")
include("init_A.jl")
include("transit_mat_entre!.jl")
include("transit_mat!.jl")
include("iterHJBentre!.jl")
include("iterHJB!.jl")
include("fdepodrift!.jl")
include("depodrift!.jl")
include("upwind!.jl")
include("Aentreassemble!.jl")
include("Aassemble!.jl")
include("Bassemble!.jl")
include("minusATassemble!.jl")
include("risky_jump_generator.jl")
include("bbuild!.jl")
include("sfbuild!.jl")
include("scebuild!.jl")
include("sdebuild!.jl")
include("cebuild!.jl")
include("findsteady.jl")
include("lab0update!.jl")
include("MIT_init.jl")
include("solveMIT!.jl")
include("solveMIT.jl")
include("MITupdate!.jl")
include("MITall!.jl")
include("initial_q_change!.jl")
include("change_distr.jl")
include("calibration_targets.jl")
include("calibrate.jl")
include("calibonlyE!.jl")
include("calibobjective.jl")
include("calibstats.jl")
include("plotMIT!.jl")
include("plotMITconv.jl")
include("shock.jl")
include("simulate_states.jl")
include("debackersim.jl")
include("makedebacker.jl")
include("percshare.jl")
include("stats.jl")
include("statsMIT.jl")
include("partial.jl")
include("checkpoint_compatibility.jl")
include("bequestregimes!.jl")
include("findsteady_fixedB.jl")
include("beqStats.jl")
include("phillips_curve.jl")
include("check_generator.jl")
include("worker_nonresponse_by_joint_quantile.jl")
include("main_fct.jl")

end
