struct TwoElectronIntegralMO{T<:Number} <: AbstractTensor
    _eri_mo::Array{T, 4}
    _fs    ::Array{T, 2}
end

function TwoElectronIntegralMO(the_scf::RestrictedSCFSolver{T,RealType}) where {T,RealType}
    if not(the_scf.is_converged)
        error("Not converged!")
    end

    nao      = get_nao(the_scf)
    nmo      = get_nmo(the_scf)
    coeff    = the_scf.orb_coeff

    temp_int = Array{T,4}(undef, nao, nao, nao, nao)
    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    @inbounds temp_int[lm, sgm, mu, nu] = get_value(
                        the_scf._fock_builder._eri, lm, sgm, mu, nu
                        )
                end
            end
        end
    end

    eri_mo = zeros(T, nao,nao,nao,nao)
    temp   = zeros(T, nao,nao,nao,nao)
    temp2  = zeros(T, nao,nao,nao,nao)
    temp3  = zeros(T, nao,nao,nao,nao)

    for mu in 1:nao
        for p in 1:nmo
            @inbounds temp[:, :, :, p] += coeff[mu, p] * temp_int[:,:,:,mu]
            for nu in 1:nao
end