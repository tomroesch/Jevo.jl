using CSV, DataFrames, Distributed, Dates


date = Dates.format(Dates.today(), "yyyy_mm_dd")

if length(ARGS) == 1
    addprocs(parse(Int64, ARGS[1]))
elseif length(ARGS) > 1
    throw(ArgumentError("Only one command line argument (cores)."))
end

@everywhere  begin
    using Jedi
    using LinearAlgebra
    using Distributions
    using SharedArrays
    emat = 2 * (ones(4, 4) - Matrix{Float64}(I, 4, 4))
    N = 1000
    f0 = 25/2N
    fl = 0.25/2N
    f = fermi_fitness(f0=f0, fl=fl)
end

rho = [0, 0.1, 0.5, 1., 2, 4]
E = SharedArray{Float64, 2}(length(rho), 1000)
L = SharedArray{Float64, 2}(length(rho), 1000)
RHO = SharedArray{Float64, 2}(length(rho), 1000)
@everywhere function run(rho, nu, l)
    pop = driver_trailer_l(N=1000, l_0=l, L=30)
    initiate_rand!(pop, 1)
    rand_rho = rand(5000000)
    rand_nu = rand(5000000)
    for i in 1:5000000
        Jedi.bp_substitution!(pop, emat, f)
        if rand_rho[i] < rho/N
            driver_mutation!(pop)
        end
        if rand_nu[i] < nu
            Jedi.l_substitution!(pop, emat, f)
        end
    end
    Gamma = get_energy(pop, emat)[1]
    l_arr = pop.l[1]
    return Gamma, l_arr
end

@sync @distributed for j in 1:1000
    for r in 1:length(rho)
        E[r, j], L[r, j] = run(rho[r], 0.01, 10)
        RHO[r, j] = rho[r]
    end
    println("Run $j done.")
end


df = DataFrame(gamma=[(E...)...], l=[(L...)...], rho=[(RHO...)...])
CSV.write(date*"_script2_results.csv", df)
