"""
    sparse_value_path(template, n)

Construct `n` sparse matrices with a common CSC sparsity pattern and independent
`nzval` arrays.  The returned matrices must only be changed through `nzval`;
structural sparse operations would change the pattern for the entire path.
"""
function sparse_value_path(template::SparseMatrixCSC{Tv,Ti}, n::Integer) where {Tv,Ti}
    n >= 0 || throw(ArgumentError("path length must be nonnegative"))

    shared_colptr = copy(template.colptr)
    shared_rowval = copy(template.rowval)
    path = Vector{SparseMatrixCSC{Tv,Ti}}(undef, n)

    for tt in eachindex(path)
        path[tt] = SparseMatrixCSC{Tv,Ti}(
            size(template, 1),
            size(template, 2),
            shared_colptr,
            shared_rowval,
            copy(template.nzval),
        )
    end

    return path
end

"""
    sparse_nz_positions(A, rows, cols)

Return the positions in `A.nzval` corresponding to `(rows[k], cols[k])`.
Throw if a requested coordinate is not already stored in the CSC pattern.  The
result permits value-only updates without using structurally mutating sparse
index assignment.
"""
function sparse_nz_positions(
    A::SparseMatrixCSC,
    rows::AbstractVector{<:Integer},
    cols::AbstractVector{<:Integer},
)
    length(rows) == length(cols) ||
        throw(DimensionMismatch("row and column index vectors must have equal length"))

    positions = Vector{Int}(undef, length(rows))
    stored_rows = rowvals(A)

    for kk in eachindex(rows, cols)
        row = rows[kk]
        col = cols[kk]
        found = 0
        for pos in nzrange(A, col)
            stored_row = stored_rows[pos]
            if stored_row == row
                found = pos
                break
            elseif stored_row > row
                break
            end
        end
        found == 0 && throw(ArgumentError(
            "coordinate ($row, $col) is not stored in the sparse pattern",
        ))
        positions[kk] = found
    end

    return positions
end
