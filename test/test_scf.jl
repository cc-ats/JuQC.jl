include("../src/JuQC.jl")
using .JuQC

int_path = "H2O_STO3G"
e_nuc, num_elec, nbas = read_data(int_path)
ovlp  = read_ovlp(int_path)
hcore = read_hcore(int_path)
eri   = read_eri(int_path)

the_scf = build_scf_solver(nbas, e_nuc, num_elec, ovlp, hcore, eri, is_restricted=true)
kernel!(the_scf, max_iter=200, tol=1e-8, the_opt=Roothaan())