module MyGlobalOptim

using Setfield, JLD2, HDF5, Sobol, Distributed, CSV, DataFrames, Optimization, OptimizationNLopt, Random, Dates, Printf, Statistics, Plots

export myglobaloptim, sobolgen, pretesting!, initmainstage, mainstage!, printinfo, resetinuse!, with_file_lock, read_locked_h5, write_locked_h5, read_locked_jld2, write_locked_jld2, update_best_if_improved!, optiminfo

"""
    myglobaloptim(objective, lbn, ubn, N, Nstar; doreset=false,
        doresetinuse=false, doresetfailed=false, dir_name="../temp/", maxtimeLocal=Inf,
        maxitersLocal=150, reltolLocal=1e-4, xtolrelLocal=1e-3,
        initial_points=nothing)

Perform a global optimization on the given `objective` function within the bounds `lbn` and `ubn`.

# Arguments
- `objective::Function`: The objective function to be minimized. If it returns a scalar, it is assumed to be the value to minimize. If it returns a tuple, it is assumed to return `(value, flag)`, where `flag` can be used to indicate if the evaluation was successful or not.
- `lbn::Vector{Float64}`: Lower bounds for the variables.
- `ubn::Vector{Float64}`: Upper bounds for the variables.
- `N::Int`: Number of initial samples.
- `Nstar::Int`: Number of samples for the main optimization stage.
- `doreset::Bool`: If `true`, doreset the optimization process and start from scratch. Default is `false`.
- `doresetinuse::Bool`: If `true`, reset the "in use" flag for Sobol points that are in use for local optimization. Default is `false`.
- `doresetfailed::Bool`: If `true`, retry failed or time-limited local searches
  without regenerating or reevaluating the Sobol design. Default is `false`.
- `dir_name::String`: Directory name for saving results. Default is `../temp/`.
- `maxtimeLocal::Int`: Maximum time for a local optimization. Default is Inf.
- `maxitersLocal::Int`: Maximum number of objective evaluations per local optimization.
- `reltolLocal::Real`: Relative objective tolerance for local optimization.
- `xtolrelLocal::Real`: Relative parameter tolerance for local optimization.
- `initial_points`: Optional vector or matrix of additional candidates, expressed
  in the optimizer's coordinates, that are appended to the Sobol design.

# Returns
- The best solution found by the optimization process.

# Description
This function performs a global optimization using the TikTak algorithm described in  Arnoud, Guvenen, and Kleineberg (2022).
1. The function first generates Sobol points in the given bounds (if `resest = false``, it reads the existing Sobol points from disk).
2. Evaluate Sobol points that have not been evaluated yet (pretesting). Results are saved to file `sobolall.h5`.
3. Select the best `Nstar` Sobol points and save them to disk as `sobolbest.h5`.
4. Perform local optimization starting from convex combination of each Sobol point with the currently best solution. File `best.jld2` is updated with the best solution found.

"""
function myglobaloptim(objective, lbn, ubn, N, Nstar; doreset=false, doresetinuse=false,
        doresetfailed=false, dir_name="../temp/", maxtimeLocal=Inf, maxitersLocal=150,
        reltolLocal=1e-4, xtolrelLocal=1e-3, initial_points=nothing)

    if !isdir(dir_name)
        mkpath(dir_name)
    end

    if doreset
        println("Resetting optimization process, creating new Sobol points and info files")
        info = (starttime=time(),
            endtime=NaN,
            finishedpre=false,
            finished=false,
            N=N,
            Nstar=Nstar,
            lbn=lbn,
            ubn=ubn,
        )

        best = (iter=0,
            x=NaN,
            val=Inf,
        )

        write_locked_jld2(dir_name * "info.jld2", info)
        write_locked_jld2(dir_name * "best.jld2", best)
        sobolgen(lbn, ubn, N, doreset, dir_name; initial_points=initial_points)
    else
        if !isfile(dir_name * "info.jld2")
            println("No info file found in " * dir_name * ", creating new one")
            info = (starttime=time(),
                endtime=NaN,
                finishedpre=false,
                finished=false,
                N=N,
                Nstar=Nstar,
                lbn=lbn,
                ubn=ubn,
            )
            write_locked_jld2(dir_name * "info.jld2", info)
        else
            info = read_locked_jld2(dir_name * "info.jld2")
        end
        if !isfile(dir_name * "best.jld2")
            println("No best file found in " * dir_name * ", creating new one")
            best = (iter=0,
                x=NaN,
                val=Inf,
            )
            write_locked_jld2(dir_name * "best.jld2", best)
        end
        sobolgen(lbn, ubn, N, doreset, dir_name; initial_points=initial_points)
        @set! info.finished = false
        @set! info.N = N
        @set! info.Nstar = Nstar
        @set! info.lbn = lbn
        @set! info.ubn = ubn
        write_locked_jld2(dir_name * "info.jld2", info)
    end
    pretesting!(objective, dir_name)
    initmainstage(Nstar, dir_name, doreset, doresetinuse, doresetfailed)
    mainstage!(objective, dir_name, maxtimeLocal, maxitersLocal,
        reltolLocal, xtolrelLocal)
    best = read_locked_jld2(dir_name * "best.jld2")

    return best.x

end

"""
    sobolgen(lbn, ubn, N, doreset, dir_name="../temp/"; initial_points=nothing)

Generate Sobol sequence points within the given bounds and save them to an HDF5 file.

# Arguments
- `lbn::Vector{Float64}`: Lower bounds for each dimension.
- `ubn::Vector{Float64}`: Upper bounds for each dimension.
- `N::Int`: Number of Sobol points to generate.
- `doreset::Bool`: If `true`, reset the Sobol points and start from scratch. Default is `false`.
- `dir_name::String`: Directory name for saving results. Default is `../temp/`.
- `initial_points`: Optional vector or matrix of points, expressed in the optimizer's
  coordinates, that are appended to the Sobol points and therefore always evaluated.

# Description
This function generates `N` Sobol sequence points within the specified lower (`lbn`) and upper (`ubn`) bounds.
The generated points are stored in a matrix `S`. If a file named `sobolall.h5` exists in the `../temp/` directory,
the function reads the existing Sobol points from the file and appends the new points to the existing ones.
Finally, the updated Sobol points are saved back to the `sobolall.h5` file.

"""
function sobolgen(lbn, ubn, N, doreset, dir_name; initial_points=nothing)
    s = SobolSeq(lbn, ubn)
    S = fill(NaN, N, length(lbn) + 2) # points, value, status flag

    for i in 1:N
        S[i, 1:end-2] = next!(s)
    end

    if initial_points !== nothing
        points = initial_points isa AbstractVector ?
                 reshape(initial_points, 1, :) : initial_points
        @assert size(points, 2) == length(lbn) "Initial points have the wrong dimension"
        @assert all((points .>= lbn') .& (points .<= ubn')) "Initial points must lie within the bounds"
        initial_rows = hcat(Matrix(points), fill(NaN, size(points, 1), 2))
        S = vcat(S, initial_rows)
    end

    filename = dir_name * "sobolall.h5"
    # If reset or no file exists, create new sobolall.h5
    if doreset || !isfile(filename)
        if !isfile(filename)
            println("No sobolall file found in " * dir_name * ", creating new one")
        end
        h5open(filename, "w") do file
            write(file, "sobolall", S)
        end
        return
    else # Load old sobolall and append new points
        S_old = h5open(filename, "r") do file
            read(file["sobolall"])
        end
        # replace rows in S where first elements are already in S_old by the old values
        # this ensures, that we do not start a local optimization from the same point twice
        D = size(S_old, 2)
        S_old_dict = Dict{NTuple{D - 2,eltype(S_old)},Vector{eltype(S_old)}}()

        for row in eachrow(S_old)
            key = Tuple(row[1:D-2])
            S_old_dict[key] = collect(row)
        end

        # Replace rows in S if key exists in S_old_dict
        for i in axes(S, 1)
            key = Tuple(S[i, 1:D-2])
            if haskey(S_old_dict, key)
                S[i, :] .= S_old_dict[key]
            end
        end
        # Check if S is any different from the old one, if not, do not save
        if S != S_old

            lockdir = filename * ".lockdir"
            with_file_lock(lockdir) do
                h5open(filename, "w") do file
                    write(file, "sobolall", S)
                end
            end
        end
    end
end

function pretesting!(objective, dir_name)

    println("Starting pretesting")

    # Find indices with NaNs in last column, i.e. those that have not been evaluated yet
    filename = dir_name * "sobolall.h5"
    S = read_locked_h5(filename, "sobolall")
    Sremain = findall(isnan, S[:, end-1])
    if isempty(Sremain)
        println("Pretesting was finished before")
        return
    else
        println(length(Sremain), " Sobol points left to evaluate")
    end
    # Sort Sobol points such that nearest neighbors are evaluated in sequence, this may speed up evaluation, because the previous transition is used as initial guess
    # As a startindex, use a random point from the list of unevaluated points. This is to avoid overlap with pretesting running in a different process
    startindex = rand(1:length(Sremain))
    Sremainsorted = nearest_neighbor_sort(S[Sremain, 1:end-2], startindex) # Think about which metric is best, it will depend on the problem
    Sremain = Sremain[Sremainsorted]
    # Parallel evaluation of corresponding Sobol points
    pmap(i -> evallist1!(i, dir_name, objective), Sremain, batch_size=1)
    println("Finished pretesting")
    info = read_locked_jld2(dir_name * "info.jld2")
    @set! info.finishedpre = true
    write_locked_jld2(dir_name * "info.jld2", info)
end

function initmainstage(Nstar, dir_name, doreset, doresetinuse, doresetfailed=false)

    filename = dir_name * "sobolall.h5"
    S = read_locked_h5(filename, "sobolall")
    CSV.write(dir_name * "sobolall.csv", DataFrame(S, :auto); writeheader=false, delim=';', decimal=',')

    # Sort S by value in second to last column and take the best x%
    N, M = size(S)
    S .= sortslices(S, by=row -> row[end-1], dims=1)
    S_best = hcat(S[1:Nstar, 1:end-1], zeros(Nstar))
    # Last column is status indicator: 0 = not used for local optimization, 1 = successfully used for local optimization, 2 = currently in use for local optimization, 3 = used unsuccessfully for local optimization (e.g. due to time limit)

    # Load old sobolbest if doreset is false
    if !doreset
        if !isfile(dir_name * "sobolbest.h5")
            println("No sobolbest file found, creating new one")
            write_locked_h5(dir_name * "sobolbest.h5", "sobolbest", S_best)
        else
            S_best_old = read_locked_h5(dir_name * "sobolbest.h5", "sobolbest")
            # replace rows in S_best where first elements are already in S_best_old by the old values
            # this ensures, that we do not start a local optimization from the same point twice
            for i = 1:Nstar
                for j = 1:length(S_best_old[:, 1])
                    if S_best_old[j, 1:end-2] == S_best[i, 1:end-2]
                        S_best[i, :] .= S_best_old[j, :]
                        break
                    end
                end
            end
            pointsleft = findall((S_best[:, end] .!= 1.0) .& (S_best[:, end] .!= 3.0))
            inuse = findall(S_best[:, end] .== 2.0)
            failed = findall(S_best[:, end] .== 3.0)
            println(length(inuse), " points are flagged as currently in use for local optimization")
            if doresetinuse
                println("Setting in use flag to 0 for these points")
                for ii in inuse
                    S_best[ii, end] = 0.0
                end
            end
            if doresetfailed
                println("Resetting ", length(failed), " failed or time-limited local optimizations")
                for ii in failed
                    S_best[ii, end] = 0.0
                end
            end
        end
    end

    # Save to sobolbest
    write_locked_h5(dir_name * "sobolbest.h5", "sobolbest", S_best)

    # Create best.h5 with the best x, corresponding y and iteration number
    if doreset & (Nstar > 0)
        best = (x=S_best[1, 1:end-2], val=S_best[1, end-1], iter=0)
        write_locked_jld2(dir_name * "best.jld2", best)
    elseif Nstar == 0
        best = (x=S[1, 1:end-2], val=S[1, end-1], iter=0)
        write_locked_jld2(dir_name * "best.jld2", best)
    else
        best = read_locked_jld2(dir_name * "best.jld2")
        if any(isnan.(best.x))
            best = (x=S_best[1, 1:end-2], val=S_best[1, end-1], iter=0)
            write_locked_jld2(dir_name * "best.jld2", best)
        end
    end
end

function mainstage!(objective, dir_name, maxtimeLocal, maxitersLocal=150,
        reltolLocal=1e-4, xtolrelLocal=1e-3)
    # Read S_best to find indices that have not been used for local optimization yet
    filename = dir_name * "sobolbest.h5"
    S_best = read_locked_h5(filename, "sobolbest") # Read the dataset
    pointsleft = findall((S_best[:, end] .!= 1.0) .& (S_best[:, end] .!= 2.0) .& (S_best[:, end] .!= 3.0))
    inuse = findall(S_best[:, end] .== 2.0)
    if isempty(pointsleft)
        println("All " * string(length(S_best[:, end])) * " Sobol points have been used for local optimization")
    else
        println(length(pointsleft) + length(inuse), " Sobol points left to be used for local optimization")
        println(length(inuse), " of these are flagged as currently in use for local optimization")
    end

    pmap(ii -> localoptim!(ii, dir_name, objective, maxtimeLocal,
        maxitersLocal, reltolLocal, xtolrelLocal), pointsleft)

    S_best = read_locked_h5(filename, "sobolbest")
    pointsleft = findall((S_best[:, end] .!= 1.0) .& (S_best[:, end] .!= 3.0))

    if isempty(pointsleft) # only update info if no more points are left
        println("Main stage is finished")
        info = read_locked_jld2(dir_name * "info.jld2")
        @set! info.endtime = time()
        @set! info.finished = true
        write_locked_jld2(dir_name * "info.jld2", info)
    end

end

function nearest_neighbor_sort(S, startindex)
    n = size(S, 1)  # Number of points
    visited = falses(n)
    sorted_indices = Vector{Int}(undef, n)

    # Start from the first point
    sorted_indices[1] = startindex
    visited[startindex] = true

    for i in 2:n
        last_index = sorted_indices[i-1]
        # Find the closest unvisited point
        dists = sum((S[[last_index], :] .- S) .^ 2, dims=2)

        dists[visited] .= Inf  # Ignore already visited points
        sorted_indices[i], _ = argmin(dists).I
        visited[sorted_indices[i]] = true
    end

    return sorted_indices
end

function evallist1!(i, dir_name, objective)
    filename = dir_name * "sobolall.h5"
    lockdir = filename * ".lockdir"
    # Read Sobol point with index i from file
    local S_i
    with_file_lock(lockdir) do
        h5open(filename, "r") do file
            S_i = file["sobolall"][i, :]
        end
    end

    if isnan(S_i[end])
        S_i = S_i[1:end-2]
    else
        println("Sobol point ", i, " has already been evaluated")
        return
    end
    println("Evaluating Sobol point ", i)
    result = objective(S_i)  # Evaluate Sobol point
    if result isa Tuple
        result, flag = result
    else
        flag = NaN
    end

    with_file_lock(lockdir) do
        h5open(filename, "r+") do file
            if isnan(file["sobolall"][i, end-1])
                file["sobolall"][i, end-1] = result
                file["sobolall"][i, end] = flag
                flush(file)
            end
        end
    end
    GC.gc()  # Force garbage collection to free memory
end

function localoptim!(ii, dir_name, objective, maxtimeLocal, maxitersLocal=150,
        reltolLocal=1e-4, xtolrelLocal=1e-3)
    # read Sobol points from row i of sobolbest
    filename = dir_name * "sobolbest.h5"
    lockdir = filename * ".lockdir"
    S_i_all = read_locked_h5_row(filename, "sobolbest", ii)

    # Check if this Sobol point has already been used for local optimization, if yes, return
    if S_i_all[end] == 1
        println("Sobol point ", ii, " has already been used for local optimization")
        return
    elseif S_i_all[end] == 2
        println("Sobol point ", ii, " currently being used for local optimization")
        return
    else # Mark as "in use" (2)
        with_file_lock(lockdir) do
            h5open(filename, "r+") do file
                file["sobolbest"][ii, end] = 2  # Separate write operation
                flush(file)
            end
        end

    end
    # The last and second to last column of S_i_all are the function value and indicator if the point has been used for local optimization, respectively
    # Consider only the Sobol point, S_i
    S_i = S_i_all[1:end-2]
    println("Local optimization starting from : ", S_i,
        " maximum evaluations: ", maxitersLocal,
        ", maximum time: ", maxtimeLocal, " seconds.")

    best = read_locked_jld2(dir_name * "best.jld2")
    info = read_locked_jld2(dir_name * "info.jld2")

    # Convex combination of Sobol point and currently best point as starting point for local optimization
    if (ii == 1) || (any(isnan.(best.x)))
        x0 = S_i
    else
        θ = min(max(0.1, sqrt(ii / info.Nstar)), 0.995)
        x0 = S_i .+ θ .* (best.x .- S_i)
    end

    # Local optimization

    f = OptimizationFunction((x, p) -> objective(x)[1])
    prob = OptimizationProblem(f, x0, 1.0, lb=info.lbn, ub=info.ubn)

    local sol
    try
        sol = Optimization.solve(prob, NLopt.LN_BOBYQA();
            maxtime=maxtimeLocal, maxiters=maxitersLocal,
            reltol=reltolLocal, xtol_rel=xtolrelLocal)
    catch err
        with_file_lock(lockdir) do
            h5open(filename, "r+") do file
                file["sobolbest"][ii, end] = 3
                flush(file)
            end
        end
        @warn "Local optimization failed" index=ii exception=(err, catch_backtrace())
        GC.gc()
        return
    end
    val = sol.objective
    x = sol.u

    println("Local optimization stopped with ", sol.retcode,
        " at ", x, " with value ", val)

    # Mark as used for local optimization in sobolbest
    completed = sol.retcode in (ReturnCode.Success, ReturnCode.MaxIters)
    with_file_lock(lockdir) do
        h5open(filename, "r+") do file
            file["sobolbest"][ii, end] = completed ? 1 : 3
            flush(file)
        end
    end

    usable = sol.retcode in (ReturnCode.Success, ReturnCode.MaxIters, ReturnCode.MaxTime)
    if usable && isfinite(val)
        best = update_best_if_improved!(dir_name, x, val, ii)
        if best.iter == ii
            println("New best solution found at  ", x, " with value ", val)
        end
    end
    GC.gc()
end

"""
resetinuse!(dir_name="../temp/")

This function resets the "in use" flag for Sobol points that are currently flagged as in use for local optimization.
"""
function resetinuse!(dir_name="../temp/")
    filename = dir_name * "sobolbest.h5"
    h5open(filename, "r+") do file
        S_best = read(file["sobolbest"])
        pointsleft = findall(S_best[:, end] .!= 1.0)
        inuse = findall(S_best[:, end] .== 2.0)
        println(length(inuse), " points are flagged as currently in use for local optimization")
        println("Setting in use flag to 0 for these points")
        for ii in inuse
            file["sobolbest"][ii, end] = 0 # Separate write operation
            flush(file)
        end
    end

end

"""
optiminfo(dir_name="../temp/")

Computes statistics about the optimization process.
"""
function optiminfo(dir_name="../temp/")
    if !isdir(dir_name)
        println("Directory " * dir_name * " does not exist")
    else
        info = read_locked_jld2(dir_name * "info.jld2")
    end
    # Pretesting info
    filename = dir_name * "sobolall.h5"
    PreTotal = NaN
    PreAdditional = NaN
    PreNoEval = NaN
    PrePercEval = NaN
    PreInf = NaN
    PrePercInf = NaN
    PreBest = NaN
    PreBestVal = NaN
    PreErrflag_counts = Dict{Float64,Int}()
    PreErrflag_shares = Dict{Float64,Float64}()
    MainUsed = NaN
    MainTimeLimit = NaN
    MainTimeLimitPerc = NaN
    MainRemaining = NaN
    MainRemainingPerc = NaN
    MainInUse = NaN
    MainInUsePerc = NaN
    MainNoEval = NaN
    MainNoEvalPerc = NaN
    try
        Sall = read_locked_h5(filename, "sobolall")
        Sall .= sortslices(Sall, by=row -> row[end-1], dims=1)
        # Write to sobolall.csv
        CSV.write(dir_name * "Allpoints.csv", DataFrame(Sall, :auto); writeheader=false, delim=';', decimal=',')
        values = Sall[:, end-1]
        evaluated = .!isnan.(values)
        PreTotal = size(Sall, 1)
        PreAdditional = max(PreTotal - info.N, 0)
        PreNoEval = sum(evaluated)
        PrePercEval = PreTotal > 0 ? PreNoEval / PreTotal * 100 : NaN
        PreInf = sum(values[evaluated] .== Inf)
        PrePercInf = PreNoEval > 0 ? PreInf / PreNoEval * 100 : NaN
        evaluated_flags = Sall[evaluated, end]
        PreErrflags = unique(evaluated_flags)
        for PreErrflag in PreErrflags
            flag_count = count(flag -> isequal(flag, PreErrflag), evaluated_flags)
            PreErrflag_counts[PreErrflag] = flag_count
            PreErrflag_shares[PreErrflag] = flag_count / PreNoEval * 100
        end
        if PreNoEval > 0
            evaluated_indices = findall(evaluated)
            best_index = evaluated_indices[argmin(values[evaluated])]
            PreBest = Sall[best_index, 1:end-2]
            PreBestVal = Sall[best_index, end-1]
        end
    catch
    end
    # Main stage info
    best = read_locked_jld2(dir_name * "best.jld2")
    if info.finishedpre
        Sbest = read_locked_h5(dir_name * "sobolbest.h5", "sobolbest")
        MainUsed = sum(Sbest[:, end] .== 1)
        MainTimeLimit = sum(Sbest[:, end] .== 3)
        MainTimeLimitPerc = sum(Sbest[:, end] .== 3) / info.Nstar * 100
        MainRemaining = sum((Sbest[:, end] .!= 1) .& (Sbest[:, end] .!= 3))
        MainRemainingPerc = MainRemaining / info.Nstar * 100
        MainInUse = sum(Sbest[:, end] .== 2)
        MainInUsePerc = sum(Sbest[:, end] .== 2) / info.Nstar * 100
        MainNoEval = info.Nstar - MainRemaining
        MainNoEvalPerc = (1.0 - MainRemaining / info.Nstar) * 100
    end

    info = (PreTotal=PreTotal,
        PreAdditional=PreAdditional,
        PreNoEval=PreNoEval,
        PrePercEval=PrePercEval,
        PreInf=PreInf,
        PrePercInf=PrePercInf,
        PreErrflag_counts=PreErrflag_counts,
        PreErrflag_shares=PreErrflag_shares,
        PreBest=PreBest,
        PreBestVal=PreBestVal,
        MainNoEval=MainNoEval,
        MainNoEvalPerc=MainNoEvalPerc,
        MainUsed=MainUsed,
        MainTimeLimit=MainTimeLimit,
        MainTimeLimitPerc=MainTimeLimitPerc,
        MainRemaining=MainRemaining,
        MainRemainingPerc=MainRemainingPerc,
        MainInUse=MainInUse,
        MainInUsePerc=MainInUsePerc,
        starttime=info.starttime,
        endtime=info.endtime,
        finishedpre=info.finishedpre,
        finished=info.finished,
        N=info.N,
        Nstar=info.Nstar,
        lbn=info.lbn,
        ubn=info.ubn,
        MainBest=best.x,
        MainBestVal=best.val
    )
    return info
end

function printinfo(dir_name="../temp/")
    # check if info.jld2 exists, if not, return
    if !isdir(dir_name)
        println("Directory " * dir_name * " does not exist")
        return
    else
        info = read_locked_jld2(dir_name * "info.jld2")
        println("###############################################")
        println("           General information")
        println("###############################################")
        println("Start time: ", Dates.format(Dates.unix2datetime(info.starttime), "yyyy-mm-dd HH:MM:SS"))
        if isnan(info.endtime)
            println("End time: Not finished")
        else
            println("End time: ", Dates.format(Dates.unix2datetime(info.endtime), "yyyy-mm-dd HH:MM:SS"))
        end
        println("Finished pretesting: ", info.finishedpre)
        println("Finished main stage: ", info.finished)
        println("Number of initial Sobol points: ", info.N)
        println("Nstar: ", info.Nstar)
        println("Lower bounds: ", info.lbn)
        println("Upper bounds: ", info.ubn)
    end

    info = optiminfo(dir_name)

    println("")
    println("###############################################")
    println("           Pretesting information")
    println("###############################################")
    println("Total pretesting points: ", info.PreTotal,
            " (", info.N, " Sobol + ", info.PreAdditional, " additional initial points)")
    if info.finishedpre
        println("Pretesting was finished before")
    else
        println("Points evaluated: ", info.PreNoEval, " (" * @sprintf("%.2f", info.PrePercEval) * "%)")

    end
    println("Points with Inf value: ", info.PreInf, " (" * @sprintf("%.2f", info.PrePercInf) * "% of evaluated points)")
    for errflag in info.PreErrflag_counts |> keys
        println("Points with error flag ", errflag, ": ", info.PreErrflag_counts[errflag], " (" * @sprintf("%.2f", info.PreErrflag_shares[errflag]) * "% of evaluated points)")
    end
    println("Best point: ", info.PreBest)
    println("Best value: ", info.PreBestVal)

    println("")
    println("###############################################")
    println("           Main stage information")
    println("###############################################")
    best = read_locked_jld2(dir_name * "best.jld2")
    println("Best solution found at: ", best.x)
    println("Best value: ", best.val)
    if info.finishedpre
        println("Sobol points used for local optimization (finished successfully): ", info.MainUsed)
        println("Sobol points used for local optimization (time limit reached): ", info.MainTimeLimit, " (", info.MainTimeLimitPerc, "%)")
        println("Sobol points remaining for local optimization (including currently running): ", info.MainRemaining, " (", info.MainRemainingPerc, "%)")
        println("Sobol points currently being used for local optimization: ", info.MainInUse, " (", info.MainInUsePerc, "%)")
    else
        println("Main stage has not started yet, still in pretesting phase")
    end
end

# The following functions are used to prevent race conditions when multiple processes try to read/write the same file.
function acquire_lock(lockpath; poll=0.05, timeout=Inf)
    t0 = time()
    while true
        try
            mkdir(lockpath)  # create directory as lock
            return nothing   # lock acquired
        catch err
            # Only treat "already exists" as "lock busy"; anything else is a real error.
            if isa(err, Base.IOError) && (err.code == -17 || err.code == Base.IOErrorCode.EEXIST)
                # If a timeout is requested, enforce it.
                if isfinite(timeout) && (time() - t0) > timeout
                    throw(ErrorException("Timeout acquiring lock $lockpath after $(timeout) seconds"))
                end

                sleep(poll)  # lock exists → wait and retry
            else
                rethrow()
            end
        end
    end
end

# Helper functions for locked HDF5 read/write
function with_file_lock(f::Function, lockpath; poll=0.05)
    acquire_lock(lockpath; poll=poll)
    try
        return f()
    finally
        rm(lockpath; force=true, recursive=true)  # release lock
    end
end

function read_locked_h5(filename::AbstractString, dataset::AbstractString)
    lockdir = filename * ".lockdir"
    return with_file_lock(lockdir) do
        h5open(filename, "r") do file
            if !haskey(file, dataset)
                error("Dataset '$(dataset)' not found in $(filename). File may be incomplete or was rewritten while another process tried to read it.")
            end
            read(file[dataset])
        end
    end
end

function read_locked_h5_row(filename::AbstractString, dataset::AbstractString, i)
    lockdir = filename * ".lockdir"
    return with_file_lock(lockdir) do
        h5open(filename, "r") do file
            if !haskey(file, dataset)
                error("Dataset '$(dataset)' not found in $(filename). File may be incomplete or was rewritten while another process tried to read it.")
            end
            file[dataset][i, :]
        end
    end
end

function write_locked_h5(filename::AbstractString, dataset::AbstractString, value)
    lockdir = filename * ".lockdir"
    with_file_lock(lockdir) do
        h5open(filename, "w") do file
            write(file, dataset, value)
            flush(file)
        end
    end
    return nothing
end

function read_locked_jld2(filename::AbstractString)
    lockdir = filename * ".lockdir"
    return with_file_lock(lockdir) do
        if !isfile(filename)
            error("File $(filename) not found.")
        end
        load_object(filename)
    end
end

function write_locked_jld2(filename::AbstractString, object)
    lockdir = filename * ".lockdir"
    with_file_lock(lockdir) do
        save_object(filename, object)
    end
    return nothing
end

function update_best_if_improved!(dir_name, x, val, iter)
    filename = dir_name * "best.jld2"
    lockdir = filename * ".lockdir"
    with_file_lock(lockdir) do
        if !isfile(filename)
            best = (x=x, val=val, iter=iter)
            save_object(filename, best)
            return best
        end
        best = load_object(filename)
        if val < best.val
            best = (x=x, val=val, iter=iter)
            save_object(filename, best)
        end
        return best
    end
end

end  # module MyGlobalOptim
