abstract type ObjectiveFunction{T} end 

function get_value(the_scf::RestrictedSCFSolver{T}) where {T}
    return the_scf.energy_elec::Real
end

function get_grad(the_scf::RestrictedSCFSolver{T}) where {T}
    return the_scf.grad::Array{T, 1}
end

abstract type OptimizationAlgorithm{T} end
function init_opt_algo!(the_opt::OptimizationAlgorithm{T}, the_obj_func::ObjectiveFunction{T}) where {T}
    
end

mutable struct RoothaanOptimizer{T} <: OptimizationAlgorithm{T}
    _is_initialized::Bool

    obj_func::Union{Nothing, ObjectiveFunction{T}}
    function RoothaanOptimizer()
        return RoothaanOptimizer{T}(false, nothing)
    end
end

function set_opt_algo!(the_opt::OptimizationAlgorithm{T}, the_scf::ObjectiveFunction{T}) where {T}
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

struct DIISOptimizer{T} <: OptimizationAlgorithm{T}
    obj_func::Union{Nothing, ObjectiveFunction{T}}
    num_sub_space::Integer
end