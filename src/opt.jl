abstract type ObjectiveFunction{T} end 

function get_value(the_scf::RestrictedSCFSolver{T})
    return the_scf.energy_elec
end

function get_grad(the_scf::RestrictedSCFSolver{T})
    return the_scf.grad
end

abstract type OptimizationAlgorithm{T} end
function init_opt_algo!(the_opt::OptimizationAlgorithm{T}, the_obj_func::ObjectiveFunction{T}) where {T}
    
end

mutable struct RoothaanOptimizer{T} <: OptimizationAlgorithm{T}
    obj_func::Union{Nothing, ObjectiveFunction{T}}
    function RoothaanOptimizer()
        return RoothaanOptimizer{T}(nothing)
    end
end

function next_step!(the_optimizer::RoothaanOptimizer{T}) where {T}
    return true
end

struct DIISOptimizer{T} <: OptimizationAlgorithm{T}
    obj_func::Union{Nothing, ObjectiveFunction{T}}
    num_sub_space::Integer
end