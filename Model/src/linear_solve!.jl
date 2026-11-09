"""
    linear_solve!(x, A, b, Temp, solver, tol; kfe=false,
                  rebuild_preconditioner=false, label="linear system")

Solve `A * x = b` with the selected linear solver. ILU factors and iterative
solver workspaces are updated and reused from `Temp`/`EntreTemp`.
"""
function linear_solve!(x, A, b, Temp, solver, tol; kfe::Bool=false,
    rebuild_preconditioner::Bool=false,
    label::AbstractString="linear system")

    if solver === :direct
        x .= A \ b
        return nothing
    end

    preconditioner = kfe ? Temp.KFEILU : Temp.HJBILU
    if rebuild_preconditioner || preconditioner === nothing
        preconditioner = ilu0(A)
        if kfe
            Temp.KFEILU = preconditioner
        else
            Temp.HJBILU = preconditioner
        end
    else
        ilu0!(preconditioner, A)
    end

    if solver === :bicgstabl
        bicgstabl!(x, A, b, 1; Pl=preconditioner, abstol=tol)
    elseif solver === :conjugategradients
        exit_code, _ = ConjugateGradients.bicgstab!(
            (out, input) -> mul!(out, A, input), b, x;
            precon=(out, input) -> ldiv!(out, preconditioner, input),
            maxIter=1000, data=Temp.BCGD, tol=tol)
        exit_code < 0 && error(
            "$label failed with :conjugategradients: " *
            ConjugateGradients.reader(exit_code))
    elseif solver === :bicgstab || solver === :dqgmres
        workspace = kfe ? Temp.KFEKrylovWork : Temp.HJBKrylovWork
        if kfe && solver === :bicgstab
            resize!(Temp.KFEKrylovShadow, length(b))
            fill!(Temp.KFEKrylovShadow, 1.0)
            Krylov.krylov_solve!(workspace, A, b, x;
                c=Temp.KFEKrylovShadow, M=preconditioner, ldiv=true,
                atol=tol, rtol=0.0)
        else
            Krylov.krylov_solve!(workspace, A, b, x;
                N=preconditioner, ldiv=true, atol=tol, rtol=0.0)
        end
        Krylov.issolved(workspace) || error(
            "$label failed with $solver: $(workspace.stats.status)")
        copyto!(x, Krylov.solution(workspace))
    else
        throw(ArgumentError(
            "Unknown solver $solver; use :bicgstabl, :bicgstab, :dqgmres, " *
            ":conjugategradients, or :direct"))
    end
    return nothing
end
