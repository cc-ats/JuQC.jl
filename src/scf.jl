abstract type FockBuilder end
struct RestrictedFockBuilder{T} <: FockBuilder
    _ovlp ::Hermitian{T,Array{T,2}}
    _hcore::Hermitian{T,Array{T,2}}
    _eri  ::TwoElectronIntegralAO{T}
end

abstract type SCFSolver{T} <: ObjectiveFunction{T} end 
mutable struct RestrictedSCFSolver{T} <: SCFSolver{T}
    _nelec::Integer
    _nocc ::Integer
    _nao  ::Integer
    _nmo  ::Integer
    _e_nuc::Real

    _fock_builder  ::RestrictedFockBuilder{T}
end

function build_scf_solver(
    nbas::Integer, e_nuc::Real, nelec::Tuple{Integer,Integer},
    s_matrix::Hermitian{T,Array{T,2}},
    h_matrix::Hermitian{T,Array{T,2}},
    eri::TwoElectronIntegralAO{T};
    is_restricted::Bool=true
    ) where {T}
    if is_restricted
        if nelec[1]==nelec[2]
            num_elec      = nelec[1] + nelec[2]
        else
            error("RHF must be closed shell")
        end
        fock_builder  = RestrictedFockBuilder{T}(s_matrix, h_matrix, eri)
        the_scf       = RestrictedSCFSolver{T}(num_elec, div(num_elec,2), nbas, nbas, e_nuc, fock_builder)
        return the_scf
    else

    end
end

function get_nao(the_scf::RestrictedSCFSolver)
    return the_scf._nao::Integer
end

function get_nmo(the_scf::RestrictedSCFSolver)
    return the_scf._nmo::Integer
end

function get_occ_index(the_scf::RestrictedSCFSolver)
    return (1:the_scf._nocc)
end

function get_vir_index(the_scf::RestrictedSCFSolver)
    return ((the_scf._nocc+1):(the_scf._nmo))
end

function get_eri_value(the_scf::RestrictedSCFSolver, lm::Integer, sgm::Integer, mu::Integer, nu::Integer)
    return get_value(the_scf._fock_builder._eri, lm, sgm, mu, nu)
end

function get_dm_tot(the_scf::RestrictedSCFSolver{T}, dm::Array{T, 2}) where {T}
    return dm
end

function get_dm_spin(the_scf::RestrictedSCFSolver{T}, dm::Array{T, 2}) where {T}
    return (dm/2, dm/2)
end

function get_e_elec(the_scf::RestrictedSCFSolver{T}, fmat::Array{T, 2}, dm::Array{T, 2}) where {T}
    hmat   = the_scf._fock_builder._hcore
    tmpmat = fmat + hmat
    return dot(dm, tmpmat)/2
end

function build_orbitals(the_scf::RestrictedSCFSolver{T}, fock_matrix::Array{T, 2}) where {T}
    ovlp_matrix = the_scf._fock_builder._ovlp
    vals, vecs  = eigen(fock_matrix, ovlp_matrix)
    return vals::Array{T,1}, vecs::Array{T,2}
end

function build_density_matrix(the_scf::RestrictedSCFSolver{T}, orb_coeff::Array{T, 2}) where {T}
    nao       = get_nao(the_scf)
    occ_index = get_occ_index(the_scf)

    co        = orb_coeff[:, occ_index]
    co_t      = transpose(co)

    dm        = Array{T}(undef, nao, nao)
    mul!(dm, co, co_t, 2.0, 0.0)
    return dm
end

function build_fock_matrix(the_scf::RestrictedSCFSolver{T}, dm::Array{T, 2}) where {T}
    nao   = get_nao(the_scf)
    hmat  = the_scf._fock_builder._hcore
    jmat  = zeros(T, nao, nao)
    kmat  = zeros(T, nao, nao)
    dm_tot            = get_dm_tot(the_scf,  dm)
    dm_alpha, dm_beta = get_dm_spin(the_scf, dm)

    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    jmat[lm, sgm] += get_eri_value(the_scf, lm, sgm, mu, nu) * dm_tot[nu, mu]
                    kmat[lm, sgm] += get_eri_value(the_scf, lm, mu, sgm, nu) * dm_alpha[nu, mu]
                end
            end
        end
    end

    fmat = hmat + jmat - kmat
    return fmat
end

function kernel(the_scf::SCFSolver{T}; max_iter::Integer=100, tol::Number=1e-8) where {T}
    
    hmat::Array{T, 2}  = the_scf._fock_builder._hcore
    orb_ene, orb_coff  = build_orbitals(the_scf, hmat)
    density_matrix     = build_density_matrix(the_scf, orb_coff)

    iter::Integer      = 0
    is_converged::Bool = false

    pre_e_elec = zero(T)
    cur_e_elec = zero(T)

    while not(is_converged) && iter < max_iter
        
        fock_matrix       = build_fock_matrix(the_scf, density_matrix)
        orb_ene, orb_coff = build_orbitals(the_scf, fock_matrix)
        
        density_matrix    = build_density_matrix(the_scf, orb_coff)
        e_elec            = get_e_elec(the_scf, fock_matrix, density_matrix)

        if not(iter == 0)
            pre_e_elec = cur_e_elec
            cur_e_elec = e_elec
            err = abs(cur_e_elec - pre_e_elec)
            is_converged = err < tol

            e_tot = e_elec  + the_scf._e_nuc
            @printf("%5d%19.10f%14.2e\n", iter, e_tot, err)
        end
        
        iter += 1
    end
end