ENV["GKSwstype"] = "100"
ENV["JULIA_WORKER_TIMEOUT"] = "600"

using Distributed
addprocs(2)

@everywhere include(joinpath(@__DIR__, "HANKModel.jl"))
@everywhere using .HANKModel

HANKModel.main_fct()
