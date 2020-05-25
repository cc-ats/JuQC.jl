using LinearAlgebra

include("../src/JuQC.jl")
using .JuQC

hm = HubbardModel(5, 100.0, 0.0)
init_guess = Matrix{Float64}(1.0*I, hm.dim, hm.dim)
this_scf = SelfConsistentField(hm, 1e-20, 100, orbtype="RHF")
c_init  = eigen(this_scf.hcore).vectors
dm_init = JuQC.get_rdm1(this_scf, c_init, this_scf.mo_occ)
kernel!(this_scf, dm_init)
