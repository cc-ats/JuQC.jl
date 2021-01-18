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

int_path = "H2O_STO3G"
e_nuc, num_elec, nbas = read_data(int_path,  FloatType=Float64, IntType=Int64)
ovlp                  = read_ovlp(int_path,  FloatType=Float64, IntType=Int64)
hcore                 = read_hcore(int_path, FloatType=Float64, IntType=Int64)
eri                   = read_eri(int_path,   FloatType=Float64, IntType=Int64)

the_scf    = build_scf_solver(nbas, e_nuc, num_elec, ovlp, hcore, eri, is_restricted=true)
scf_result = kernel!(the_scf, max_iter=200, tol=1e-8, scf_algo=DIIS(T=Float64))
get_mp2_ecorr(scf_result)

