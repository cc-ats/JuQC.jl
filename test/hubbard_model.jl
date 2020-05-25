include("../src/JuQC.jl")
using .JuQC

hm = HubbardModel(5, 1.0, 4.0)
init_guess = Matrix{Float64}(1.0*I, hm.dim, hm.dim)
this_scf = SelfConsistentField(hm, 1e-20, 10, init_guess, orbtype="UHF")
kernel!(this_scf)