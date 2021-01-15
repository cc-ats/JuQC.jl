abstract type SCFAlgorithm{T} end

mutable struct Roothaan{T} <: SCFAlgorithm{T}
    _is_initialized::Bool
    the_scf::Union{Nothing, SCFSolver{T}}
end

function Roothaan(;T=Float64)
    return Roothaan{T}(false, nothing)
end

function init_scf_algo!(the_opt::SCFAlgorithm{T}, the_scf::ObjectiveFunction{T}) where {T}
    if isnothing(the_opt.the_scf) && not(the_opt._is_initialized)
        the_opt._is_initialized = true
        the_opt.the_scf        = the_scf
    else
        error("The objective function is not nothing.")
    end
end

function next_step!(the_optimizer::Roothaan{T}) where {T}
    return true
end

function get_density_matrix(the_optimizer::Roothaan{T}) where {T}
    if not(isnothing(the_optimizer.the_scf))
        return the_optimizer.the_scf.density_matrix
    end
end

struct DIIS{T} <: SCFAlgorithm{T}
    _is_initialized::Bool
    _num_sub_space::Integer
    the_scf::Union{Nothing, ObjectiveFunction{T}}
end

function DIIS(;num_sub_space=8, T=Float64)
    return DIIS{T}(false, num_sub_space, nothing)
end

function init_scf_algo!(the_opt::SCFAlgorithm{T}, the_scf::ObjectiveFunction{T}) where {T}
    if isnothing(the_opt.the_scf) && not(the_opt._is_initialized)
        the_opt._is_initialized = true
        the_opt.the_scf        = the_scf
    else
        error("The objective function is not nothing.")
    end
end

function next_step!(the_optimizer::Roothaan{T}) where {T}
    return true
end

function get_density_matrix(the_optimizer::Roothaan{T}) where {T}
    if not(isnothing(the_optimizer.the_scf))
        return the_optimizer.the_scf.density_matrix
    end
end