struct TwoElectronIntegralMO{T<:Number} <: AbstractTensor
    _eri_mo_spin::Array{T, 4}
    _f_mo_spin  ::Array{T, 1}
end

function build_eri_mo(the_scf::RestrictedSCFResult{T,RealType}) where {T,RealType}
    if not(the_scf.is_converged)
        error("Not converged!")
    end

    nao      = get_nao(the_scf)
    nmo      = get_nmo(the_scf)
    coeff    = get_orb_coeff(the_scf)
    ene      = get_orb_ene(the_scf)

    temp_int = Array{T,4}(undef, nao, nao, nao, nao)
    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    @inbounds temp_int[lm, sgm, mu, nu] = get_value(
                        the_scf.fock_builder._eri, lm, sgm, mu, nu
                        )
                end
            end
        end
    end

    t0 = time()
    
    temp1  = zeros(T, nao,nao,nao,nmo)
    temp2  = zeros(T, nao,nao,nmo,nmo)
    temp3  = zeros(T, nao,nmo,nmo,nmo)
    eri_mo = zeros(T, nmo,nmo,nmo,nmo)

    for p in 1:nmo
        for mu in 1:nao
            @inbounds temp1[p,:,:,:] += coeff[mu, p] * temp_int[mu,:,:,:]
        end
        for q in 1:nmo
            for nu in 1:nao
                @inbounds temp2[p,q,:,:] += coeff[nu, q] * temp1[p,nu,:,:]
            end
            for r in 1:nmo
                for lm in 1:nao
                    @inbounds temp3[p,q,r,:] += coeff[lm, r] * temp2[p,q,lm,:]
                end
                for s in 1:nmo
                    for sg in 1:nao
                        @inbounds eri_mo[p,q,r,s] += coeff[sg, s] * temp3[p,q,r,sg]
                    end
                end
            end
        end
    end

    eri_mo_spin = Array{T,4}(undef, 2*nmo, 2*nmo, 2*nmo, 2*nmo)

    for p in 1:2*nmo
        for q in 1:2*nmo
            for r in 1:2*nmo
                for s in 1:2*nmo
                    pp = div((p-1), 2) + 1
                    qq = div((q-1), 2) + 1
                    rr = div((r-1), 2) + 1
                    ss = div((s-1), 2) + 1

                    @inbounds val1::T = eri_mo[pp, rr, qq, ss] * (p%2 == r%2) * (q%2 == s%2)
                    @inbounds val2::T = eri_mo[pp, ss, qq, rr] * (p%2 == s%2) * (q%2 == r%2)
                    @inbounds eri_mo_spin[p, q, r, s] = (val1 - val2)::T

                end
            end
        end
    end

    f_spin = zeros(T, 2*nmo)
    for p in 1:2*nmo
        f_spin[p] = ene[div((p-1),2)+1]::T
    end

    println("Time for AO2MO is ", time()-t0, "\n")

    return TwoElectronIntegralMO{T}(eri_mo_spin, f_spin)
end

function get_mp2_ecorr(the_scf::RestrictedSCFResult{T,RealType}) where {T,RealType}
    eri_mo = build_eri_mo(the_scf)

    occ = get_occ_index(the_scf)
    vir = get_vir_index(the_scf)

    e_corr      = zero(RealType)
    eri_mo_spin = eri_mo._eri_mo_spin
    f_mo_spin   = eri_mo._f_mo_spin

    for i in occ
        for j in occ
            for a in vir
                for b in vir
                    ia = 2*(i-1)+1
                    ib = 2*(i-1)+2

                    ja = 2*(j-1)+1
                    jb = 2*(j-1)+2

                    aa = 2*(a-1)+1
                    ab = 2*(a-1)+2

                    ba = 2*(b-1)+1
                    bb = 2*(b-1)+2

                    e_corr += 0.25*(eri_mo_spin[ia,ja,aa,ba] * eri_mo_spin[ia,ja,aa,ba])/(f_mo_spin[ia] + f_mo_spin[ja] - f_mo_spin[aa] - f_mo_spin[ba])

                    e_corr += 0.25*(eri_mo_spin[ia,jb,ab,ba] * eri_mo_spin[ia,jb,ab,ba])/(f_mo_spin[ia] + f_mo_spin[jb] - f_mo_spin[ab] - f_mo_spin[ba])

                    e_corr += 0.25*(eri_mo_spin[ib,ja,aa,bb] * eri_mo_spin[ib,ja,aa,bb])/(f_mo_spin[ib] + f_mo_spin[ja] - f_mo_spin[aa] - f_mo_spin[bb])

                    e_corr += 0.25*(eri_mo_spin[ia,jb,aa,bb] * eri_mo_spin[ia,jb,aa,bb])/(f_mo_spin[ia] + f_mo_spin[jb] - f_mo_spin[aa] - f_mo_spin[bb])

                    e_corr += 0.25*(eri_mo_spin[ib,ja,ab,ba] * eri_mo_spin[ib,ja,ab,ba])/(f_mo_spin[ib] + f_mo_spin[ja] - f_mo_spin[ab] - f_mo_spin[ba])

                    e_corr += 0.25*(eri_mo_spin[ib,jb,ab,bb] * eri_mo_spin[ib,jb,ab,bb])/(f_mo_spin[ib] + f_mo_spin[jb] - f_mo_spin[ab] - f_mo_spin[bb])
                    
                end
            end
        end
    end
    return e_corr::RealType
end
