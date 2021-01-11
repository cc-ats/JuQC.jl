abstract type ObjectiveFunction{T} end 

abstract type OptimizationAlgorithm{T} end

mutable struct RoothaanOptimizer{T} <: OptimizationAlgorithm{T}
    _the_obj_func::ObjectiveFunction{T}

    _cur_y::Real
    _pre_y::Real

    _cur_grad::Array{T,1}
    _pre_grad::Array{T,1}

    _cur_x::Array{T,1}
    _pre_x::Array{T,1}

    _the_scf::RestrictedSCFSolver{T}
end

function next_step!(the_optimizer::RoothaanOptimizer{T}, this_y::T) where {T}
    the_optimizer._prev_y    = the_optimizer._this_y
    the_optimizer._this_y    = this_y

    grad_size                = sizeof(the_optimizer._grad_norm)
    the_optimizer._grad_norm = norm(get_grad(the_optimizer._the_obj_func))/sqrt(grad_size)
    
    is_wolfe1 = abs(_this_y-_prev_y)     < the_optimizer._tol
    is_wolfe2 = the_optimizer._grad_norm < the_optimizer._tol

    return is_wolfe1 && is_wolfe2
end

struct DIISOptimizer{T}     <: OptimizationAlgorithm{T}
    tol::Number
end