### This is the main file for the model part of "Monetary Policy and Aggregate
### Investment: The Role of Entrepreneurs" by A. Matusche & J.
### Wacks

ENV["GKSwstype"] = "100"

nothreads = Threads.nthreads()
println("Number of threads: " * string(nothreads))

include("HANKModel.jl")
using .HANKModel

HANKModel.main_fct()
