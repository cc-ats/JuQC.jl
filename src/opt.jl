abstract type ObjectiveFunction{T} end

abstract type OptimizationAlgorithm{T} end

mutable struct RoothaanOptimizer{T} <: OptimizationAlgorithm{T}
    _is_initialized::Bool
    obj_func::Union{Nothing, ObjectiveFunction{T}}
end

function RoothaanOptimizer(;T=Float64)
    return RoothaanOptimizer{T}(false, nothing)
end

function init_opt_algo!(the_opt::OptimizationAlgorithm{T}, the_scf::ObjectiveFunction{T}) where {T}
    if isnothing(the_opt.obj_func) && not(the_opt._is_initialized)
        the_opt._is_initialized = true
        the_opt.obj_func        = the_scf
    else
        error("The objective function is not nothing.")
    end
end

function next_step!(the_optimizer::RoothaanOptimizer{T}) where {T}
    return true
end

function get_density_matrix(the_optimizer::RoothaanOptimizer{T}) where {T}
    if not(isnothing(the_optimizer.obj_func))
        return the_optimizer.obj_func.density_matrix
    end
end

struct DIISOptimizer{T} <: OptimizationAlgorithm{T}
    obj_func::Union{Nothing, ObjectiveFunction{T}}
    num_sub_space::Integer
end