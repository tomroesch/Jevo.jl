using Distributed
#addprocs(48)

@everywhere  begin
    using Jedi
    using LinearAlgebra
    using Distributions
    using SharedArrays
    emat = 2 * (ones(4, 4) - Matrix{Float64}(I, 4, 4))
    Est(l) = 2*(3/4 * l - 5)
    N = 1000
    f0 = 20/2N
    f = fermi_fitness(f0=f0, E_Star=Est)
end

results = SharedArray{Float64, 3}(2, 4, 100)

rho = [0, 0.1, 0.5, 1., 2]
E = SharedArray{Float64, 2}(length(rho), 1000)
L = SharedArray{Float64, 2}(length(rho), 1000)
RHO = SharedArray{Float64, 2}(length(rho), 1000)



@everywhere function run(rho, nu, l)
    pop = driver_trailer_l(N=1000, l_0=l, L=l)
    initiate_rand!(pop, 1)
    for i in 1:5000000
        Jedi.bp_substitution!(pop, emat, f)
        if rand() < rho/N
            driver_mutation!(pop)
        end
        if rand() < nu/N
            Jedi.l_substitution!(pop, emat, f)
        end
    end
    Gamma = sum(get_energy(pop, emat) .* pop.freqs ./ pop.l) / pop.N
    l_arr = sum(pop.l .* pop.freqs) / pop.N
    return Gamma, l_arr
end

@sync @distributed for j in 1:1000
    for r in 1:length(rho)
        E[r, j], L[r, j] = run(rho[r], 0.01, l)
        L[r, j] = pop.l[1]
        RHO[r, j] = rho[r]
    end
    println("Run $j done.")
end


df = DataFrame(gamma=[(E...)...], l=[(L...)...], rho=[(RHO...)...])
CSV.write("script2_results.csv", df)
