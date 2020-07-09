"""
    function initiate_rand!(pop::binding_sites, c)

Iniate a population of c species of random sequences.
"""
function initiate_rand!(pop::binding_sites, c)
    N_sub = pop.N ÷ c
    rest = pop.N - N_sub * c
    for i in 1:c
        push!(pop.seqs, rand(collect(1:pop.n), pop.L))
        push!(pop.frequs, N_sub)
    end
    pop.frequs[1] += rest
end
    