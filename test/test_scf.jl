include("../src/JuQC.jl")
using .JuQC

int_path = "H2O_STO3G"
e_nuc, num_elec, nbas = read_data(int_path,  FloatType=Float64, IntType=Int64)
ovlp                  = read_ovlp(int_path,  FloatType=Float64, IntType=Int64)
hcore                 = read_hcore(int_path, FloatType=Float64, IntType=Int64)
eri                   = read_eri(int_path,   FloatType=Float64, IntType=Int64)

the_scf    = build_scf_solver(nbas, e_nuc, num_elec, ovlp, hcore, eri, is_restricted=true)
scf_result = kernel!(the_scf, max_iter=200, tol=1e-8, scf_algo=DIIS(T=Float64))

println("\n#########################################")
println("scf_converged = ", scf_result.is_converged)
println("energy_tot  = ", scf_result.energy_tot)
println("energy_elec = ", scf_result.energy_elec)
println("#########################################")

println("\n\norb_ene = ")
display(scf_result.orb_ene)

println("\n\norb_coeff = ")
display(scf_result.orb_coeff)

println("\n\ngrad = ")
display(scf_result.grad)

println("\n\ndensity_matrix = ")
display(scf_result.density_matrix)

println("\n\nfock_matrix = ")
display(scf_result.fock_matrix)
