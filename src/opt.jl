abstract type ObjectiveFunction{T} end 

abstract type OptimizationAlgorithm{T} end

mutable struct RoothaanOptimizer{T} <: OptimizationAlgorithm{T}
    _the_obj_func::ObjectiveFunction{T}

    _tol   ::T
    _prev_y::T
    _this_y::T
    _grad_norm::T
end

function next_step!(the_optimizer::RoothaanOptimizer{T}, this_y::T) where {T}
    the_optimizer._prev_y    = the_optimizer._this_y
    the_optimizer._this_y    = this_y
    the_optimizer._grad_norm = norm(get_grad(the_optimizer._the_obj_func))
    
    is_wolfe1 = abs(_this_y-_prev_y)     < the_optimizer._tol
    is_wolfe2 = the_optimizer._grad_norm < the_optimizer._tol

    return is_wolfe1 && is_wolfe2
end

struct DIISOptimizer{T}     <: OptimizationAlgorithm{T}
    tol::Number
end