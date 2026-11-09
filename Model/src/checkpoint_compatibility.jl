# Solver settings describe how a path was computed, not which equilibrium it
# solves. Keep them in saved signatures for provenance, but exclude them when
# deciding whether a checkpoint can initialize another solve.
function compatible_checkpoint_signature(saved, current)
    saved isa NamedTuple && current isa NamedTuple || return false
    economic_signature(signature) = (; (key => value for (key, value) in
        pairs(signature) if key != :solver_options)...)
    return isequal(economic_signature(saved), economic_signature(current))
end

function checkpoint_solver_settings_match(saved, current)
    return saved isa NamedTuple && current isa NamedTuple &&
           hasproperty(saved, :solver_options) &&
           hasproperty(current, :solver_options) &&
           isequal(saved.solver_options, current.solver_options)
end
