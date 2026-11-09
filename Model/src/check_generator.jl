function check_generator(A, dims; tol=1e-8)
    rowsum = vec(sum(A, dims=2))
    I, J, V = findnz(A)
    offdiag = [v for (i, j, v) in zip(I, J, V) if i != j]

    max_rowsum = maximum(abs.(rowsum))
    min_offdiag = minimum(offdiag)
    max_diag = maximum(diag(A))

    if max_rowsum > tol || min_offdiag < -tol || max_diag > tol
        badrow = argmax(abs.(rowsum))
        badstate = CartesianIndices(dims)[badrow]
        @show max_rowsum min_offdiag max_diag badrow badstate
        error("Matrix is not a valid generator")
    end
end
