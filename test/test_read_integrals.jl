using Test

int_path = "H2O_STO3G"
e_nuc, num_elec, nbas = read_data(int_path)
println("e_nuc    = ", e_nuc)
println("num_elec = ", num_elec)
println("nbas     = ", nbas)

eri = read_eri(int_path)
for lm in 1:nbas
    for mu in 1:nbas
        for sgm in 1:nbas
            for rho in 1:nbas
                @test get_value(eri, lm, mu, sgm, rho) == get_value(eri, sgm, rho, lm, mu)
                @test get_value(eri, lm, mu, sgm, rho) == get_value(eri, mu, lm, sgm, rho)
                @test get_value(eri, lm, mu, sgm, rho) == get_value(eri, lm, mu, rho, sgm)
            end
        end
    end
end