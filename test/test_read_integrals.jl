include("../src/JuQC.jl")
using .JuQC

int_path = "H2O_STO3G"
e_nuc, num_elec, nbas = read_data(int_path)
println("e_nuc    = ", e_nuc)
println("num_elec = ", num_elec)
println("nbas     = ", nbas)

ovlp = read_ovlp(int_path)
println("ovlp = ")
for lm in 1:6
    for mu in 1:6
        println("get_value(ovlp, $lm, $mu) = ", get_value(ovlp, lm, mu))
    end
end

hcore = read_hcore(int_path)
println("hcore = ")
display(hcore._data)
for lm in 1:nbas
    for mu in 1:nbas
        println("get_value(hcore, $lm, $mu) = ", get_value(hcore, lm, mu))
    end
end

eri = read_eri(int_path)
println("eri = ")
for lm in 1:nbas
    for mu in 1:nbas
        for sgm in 1:nbas
            for rho in 1:nbas
                println("get_value(eri, $lm, $mu, $sgm, $rho) = ", get_value(eri, lm, mu, sgm, rho))
            end
        end
    end
end