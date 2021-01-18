struct FockBuilder{T}
    _nao  ::Integer
    _ovlp ::Hermitian{T,Array{T,2}}
    _hcore::Hermitian{T,Array{T,2}}
    _eri  ::TwoElectronIntegralAO{T}
end

abstract type  SCFResult{T,RealType<:Real} end 
struct RestrictedSCFResult{T,RealType} <: SCFResult{T,RealType}
    nao           ::Integer
    nmo           ::Integer
    nocc          ::Integer
    e_nuc         ::RealType
    fock_builder  ::FockBuilder{T}

    is_converged  ::Bool

    # results
    energy_tot    ::RealType
    energy_elec   ::RealType

    orb_ene       ::Array{RealType,1}
    orb_coeff     ::Array{T,2}
    grad          ::Array{T,2}

    density_matrix::Hermitian{T,Array{T,2}}
    fock_matrix   ::Hermitian{T,Array{T,2}}
end

function get_nao(the_scf::RestrictedSCFResult)
    return the_scf.nao::Integer
end

function get_nmo(the_scf::RestrictedSCFResult)
    return the_scf.nmo::Integer
end

function get_nocc(the_scf::RestrictedSCFResult)
    return the_scf.nocc
end

function get_nvir(the_scf::RestrictedSCFResult)
    return the_scf.nmo - the_scf.nocc
end

function get_occ_index(the_scf::RestrictedSCFResult)
    return (1:the_scf.nocc)
end

function get_vir_index(the_scf::RestrictedSCFResult)
    return ((the_scf.nocc+1):(the_scf.nmo))
end

function get_orb_ene(the_scf::RestrictedSCFResult{T,RealType}) where {T,RealType}
    return the_scf.orb_ene::Array{RealType,1}
end

function get_orb_coeff(the_scf::RestrictedSCFResult{T,RealType}) where {T,RealType}
    return the_scf.orb_coeff::Array{T,2}
end

struct UnrestrictedSCFResult{T,RealType}      <: SCFResult{T,RealType}
    is_converged  ::Bool

    # results
    energy_tot    ::RealType
    energy_elec   ::RealType

    orb_ene       ::Tuple{Array{RealType,1},Array{RealType,1}} 
    orb_coeff     ::Tuple{Array{T,2},Array{T,2}}
    grad          ::Tuple{Array{T,2},Array{T,2}}

    density_matrix::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
    fock_matrix   ::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
end

abstract type  SCFSolver{T,RealType<:Real} end 
mutable struct RestrictedSCFSolver{T,RealType} <: SCFSolver{T,RealType}
    # input parameters
    _nao           ::Integer
    _nmo           ::Integer
    _nocc          ::Integer
    _e_nuc         ::RealType
    _fock_builder  ::FockBuilder{T}

    # parameters
    is_init_guess  ::Bool
    is_converged   ::Bool

    # results
    num_fock_build ::Integer
    time_fock_build::RealType

    energy_tot     ::RealType
    energy_elec    ::RealType

    orb_ene        ::Array{RealType,1}
    orb_coeff      ::Array{T,2}
    grad           ::Array{T,2}

    density_matrix ::Hermitian{T,Array{T,2}}
    fock_matrix    ::Hermitian{T,Array{T,2}}
end

mutable struct UnrestrictedSCFSolver{T, RealType} <: SCFSolver{T,RealType}
    _nao           ::Integer
    _nmo           ::Integer
    _nocc          ::Tuple{Integer,Integer}
    _e_nuc         ::RealType
    _fock_builder  ::FockBuilder{T}

    is_init_guess  ::Bool
    is_converged   ::Bool

    num_fock_build ::Integer
    time_fock_build::RealType

    energy_tot     ::RealType
    energy_elec    ::RealType

    orb_ene        ::Tuple{Array{RealType,1},Array{RealType,1}} 
    orb_coeff      ::Tuple{Array{T,2},Array{T,2}}
    grad           ::Tuple{Array{T,2},Array{T,2}}

    density_matrix::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
    fock_matrix   ::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
end

function build_scf_solver(
    nao::Integer, e_nuc::Real,
    nelec::Tuple{Integer,Integer},
    ovlp::Hermitian{T,Array{T,2}},
    hcore::Hermitian{T,Array{T,2}},
    eri::TwoElectronIntegralAO{T};
    is_restricted::Bool=true
    ) where {T}

    if is_restricted
        if not(nelec[1]==nelec[2])
            error("RHF must be closed shell")
        end
        fock_builder  = FockBuilder{T}(nao, ovlp, hcore, eri)
        the_scf       = RestrictedSCFSolver{T,real(T)}(
            nao, nao, nelec[1], e_nuc, fock_builder,
            false, false, 0, 0.0, 0.0, 0.0,
            zeros(T, nao), zeros(T, nao, nao), zeros(T, nao, nao),
            Hermitian(zeros(T, nao, nao)), Hermitian(zeros(T, nao, nao))
            )
        return the_scf
    else
        error("Not implemented!")
    end
end

function get_nao(the_scf::RestrictedSCFSolver)
    return the_scf._nao::Integer
end

function get_nmo(the_scf::RestrictedSCFSolver)
    return the_scf._nmo::Integer
end

function get_nocc(the_scf::RestrictedSCFSolver)
    return the_scf._nocc
end

function get_nvir(the_scf::RestrictedSCFSolver)
    return the_scf._nmo - the_scf._nocc
end

function get_occ_index(the_scf::RestrictedSCFSolver)
    return (1:the_scf._nocc)
end

function get_vir_index(the_scf::RestrictedSCFSolver)
    return ((the_scf._nocc+1):(the_scf._nmo))
end

function get_density_matrix_tot(the_scf::RestrictedSCFSolver{T,RealType}, p::Hermitian{T,Array{T,2}}) where {T,RealType}
    return p::Hermitian{T,Array{T,2}}
end

function get_density_matrix_spin(the_scf::RestrictedSCFSolver{T,RealType}, p::Hermitian{T,Array{T,2}}) where {T,RealType}
    return ((p/2, p/2))::Tuple{Hermitian{T,Array{T,2}}, Hermitian{T,Array{T,2}}}
end

function get_e_elec(the_scf::RestrictedSCFSolver{T,RealType}, f::Hermitian{T,Array{T,2}}, p::Hermitian{T,Array{T,2}}) where {T,RealType}
    h      = the_scf._fock_builder._hcore
    return real(dot(p, f + h)/2)::RealType
end

function get_grad(the_scf::RestrictedSCFSolver{T,RealType}, f::Hermitian{T,Array{T,2}}, p::Hermitian{T,Array{T,2}}) where {T,RealType}
    s    = the_scf._fock_builder._ovlp
    fps  = f * p * s
    g    = fps - transpose(fps)
    n    = norm(g)
    return n::RealType, g::Array{T,2}
end

function build_orbitals(
    the_scf::RestrictedSCFSolver{T,RealType},
    fock::Hermitian{T,Array{T,2}}
    ) where {T,RealType}
    ovlp        = the_scf._fock_builder._ovlp
    vals, vecs  = eigen(
        fock::Hermitian{T,Array{T,2}},
        ovlp::Hermitian{T,Array{T,2}}
        )
    return vals::Array{RealType,1}, vecs::Array{T,2}
end

function build_density_matrix(
    the_scf::RestrictedSCFSolver{T,RealType},
    orb_coeff::Array{T, 2}
    ) where {T,RealType}
    nao       = get_nao(the_scf)
    occ_index = get_occ_index(the_scf)

    co        = orb_coeff[:, occ_index]
    co_adj    = adjoint(co)
    dm        = 2 * co * co_adj

    return Hermitian(dm::Array{T,2})
end

function build_fock_matrix(
    the_scf::RestrictedSCFSolver{T,RealType},
    dm::Hermitian{T,Array{T,2}}
    ) where {T,RealType}
    nao    = the_scf._nao
    hcore  = the_scf._fock_builder._hcore

    dm_tot            = get_density_matrix_tot(the_scf::RestrictedSCFSolver{T,RealType},  dm)
    dm_alpha, dm_beta = get_density_matrix_spin(the_scf::RestrictedSCFSolver{T,RealType}, dm)

    t0 = time()
    jmat   = zeros(T, nao, nao)
    kmat   = zeros(T, nao, nao)

    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    @inbounds jmat[lm, sgm] += get_value(
                        the_scf._fock_builder._eri, lm, sgm, mu, nu
                        ) * dm_tot[nu, mu]

                    @inbounds kmat[lm, sgm] += get_value(
                        the_scf._fock_builder._eri, lm, mu, sgm, nu
                        ) * dm_alpha[nu, mu]
                end
            end
        end
    end

    fock = hcore + jmat - kmat
    the_scf.num_fock_build  += 1
    the_scf.time_fock_build += time() - t0

    return Hermitian(fock::Array{T,2})
end

function set_init_guess!(
    the_scf::RestrictedSCFSolver{T,RealType}, 
    tmp::Array{T,2};
    is_orb_coeff::Bool=true, is_density_matrix::Bool=false
    ) where {T,RealType}
    nao = get_nao(the_scf)
    nmo = get_nmo(the_scf)

    if is_orb_coeff && not(is_density_matrix)
        if not(sizeof(tmp) == (nao, nmo))
            error("Wrong init guess dimension!!!")
        end
        the_scf.is_init_guess  = true
        the_scf.density_matrix = build_density_matrix(the_scf, tmp)
    elseif not(is_orb_coeff) && is_density_matrix
        if not(sizeof(tmp) == (nao, nao))
            error("Wrong init guess dimension!!!")
        end
        the_scf.is_init_guess  = true
        the_scf.density_matrix = Hermitian(tmp)
    else
        error("Unknown init guess type!!!")
    end
end

function update!(
    the_scf       ::RestrictedSCFSolver{T,RealType},
    is_converged  ::Bool,
    energy_tot    ::RealType,
    energy_elec   ::RealType,
    orb_ene       ::Array{T,1},
    orb_coeff     ::Array{T,2},
    grad          ::Array{T,2},
    density_matrix::Hermitian{T,Array{T,2}},
    fock_matrix   ::Hermitian{T,Array{T,2}},
    ) where {T, RealType}

    the_scf.is_converged   = is_converged
    the_scf.energy_tot     = energy_tot
    the_scf.energy_elec    = energy_elec
    the_scf.orb_ene        = orb_ene
    the_scf.orb_coeff       = orb_coeff
    the_scf.grad           = grad
    the_scf.density_matrix = density_matrix
    the_scf.fock_matrix    = fock_matrix
end