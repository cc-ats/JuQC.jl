abstract type FockBuilder end

struct RestrictedFockBuilder{T} <: ObjectiveFunction
    _ovlp ::Array{T, 2}
    _hcore::Array{T, 2}
    _eri  ::TwoElectronIntegralAO{T}
end

abstract type SCFSolver end 
mutable struct RestrictedSCFSolver{T} <: SCFSolver
    _nelec::Integer
    _nocc ::Integer
    _nao  ::Integer
    _nmo  ::Integer

    _orb_ene       ::Array{T, 1}
    _orb_coeff     ::Array{T, 2}
    _orb_coeff_inv ::Array{T, 2}
    _density_matrix::Array{T, 2}

    _guess_method  ::InitGuessMethod
    _fock_builder  ::RestrictedFockBuilder{T}
end

function get_nao(the_scf::RestrictedSCFSolver)
    return the_scf._nao::Integer
end

function get_nmo(the_scf::RestrictedSCFSolver)
    return the_scf._nmo::Integer
end

function get_orb_coeff_value(the_scf::RestrictedRestrictedSCFSolver{T}, orb_coeff::Array{T, 2}, lm::Integer, p::Integer) where {T}
    return orb_coeff[lm, p]::T
end

function get_orb_coeff_vector(the_scf::RestrictedRestrictedSCFSolver{T}, orb_coeff::Array{T, 2}, p::Integer) where {T}
    return orb_coeff[:, p]::Array{T, 1}
end

function get_orb_coeff_block(the_scf::RestrictedRestrictedSCFSolver{T}, p::UnitRange{Integer}) where {T}
    return the_scf._orb_coeff[:, p]::Array{T, 2}
end

function get_orb_ene(the_scf::RestrictedRestrictedSCFSolver{T}, p::Integer) where {T}
    return the_scf._orb_ene[p]::T
end

function get_occ_index(the_scf::RestrictedRestrictedSCFSolver)
    return (1:the_scf._nocc)::UnitRange{Integer}
end

function get_vir_index(the_scf::RestrictedRestrictedSCFSolver)
    return ((the_scf._nocc+1):(the_scf._nmo))::UnitRange{Integer}
end

function build_density_matrix(the_scf::RestrictedRestrictedSCFSolver{T}, orb_coeff::Array{T, 2}) where {T}
    co = get
end

function get_fock(the_scf::RestrictedRestrictedSCFSolver{T}, dm::Array{T, 2}) where {T}
    
    jmat = Array{T, 2}(undef, )
end

struct UnrestrictedRestrictedSCFSolver{T}
    _nelec::Tuple{Integer,Integer}
    _nocc ::Tuple{Integer,Integer}
    _nao  ::Integer
    _nmo  ::Integer

    _orb_ene       ::Array{T, 2}
    _orb_coeff     ::Array{T, 3}
    _orb_coeff_inv ::Array{T, 3}
    _density_matrix::Array{T, 3}

    _guess_method  ::InitGuessMethod
    _fock_builder  ::UnrestrictedFockBuilder{T}
end

function get_nao(the_scf::UnrestrictedSCFSolver)
    return the_scf._nao::Integer
end

function get_nmo(the_scf::UnrestrictedSCFSolver)
    return the_scf._nmo::Integer
end

function get_orb_coeff_value(the_scf::UnrestrictedUnrestrictedSCFSolver{T}, s::Integer, lm::Integer, p::Integer) where {T}
    return the_scf._orb_coeff[s, lm, p]::T
end

function get_orb_coeff_vector(the_scf::UnrestrictedUnrestrictedSCFSolver{T}, s::Integer, p::Integer) where {T}
    return the_scf._orb_coeff[s, :, p]::Array{T, 1}
end

function get_orb_coeff_block(the_scf::UnrestrictedUnrestrictedSCFSolver{T}, s::Integer, p::UnitRange{Integer}) where {T}
    return the_scf._orb_coeff[s, :, p]::Array{T, 2}
end

function get_orb_ene(the_scf::UnrestrictedUnrestrictedSCFSolver{T}, s::Integer, p::Integer) where {T}
    return the_scf._orb_ene[s, p]::T
end

function get_occ_index(the_scf::UnrestrictedUnrestrictedSCFSolver)
    return ((1:the_scf._nocc[1]), (1:the_scf._nocc[2]))::Tuple{UnitRange{Integer},UnitRange{Integer}}
end

function get_vir_index(the_scf::UnrestrictedUnrestrictedSCFSolver)
    return ((the_scf._nocc[1]+1):(the_scf._nmo), (the_scf._nocc[2]+1):(the_scf._nmo))::Tuple{UnitRange{Integer},UnitRange{Integer}}
end

function kernel()
    
end