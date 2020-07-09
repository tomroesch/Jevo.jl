using StatsBase

function mutation!(pop::binding_sites)
    mutation_sequence = sample(collect(1:length(pop.seqs)), Weights(pop.freqs ./ pop.N))
    push!(pop.seqs, copy(pop.seqs[mutation_sequence]))
    push!(pop.freqs, 1)
    pop.freqs[mutation_sequence] -= 1
    mutation_base = rand(1:pop.L)
    poss_bases = filter(x->x != pop.seqs[end][mutation_base], collect(1:pop.n))
    pop.seqs[end][mutation_base] = rand(poss_bases)
    remove_empty!(pop.freqs, pop.seqs)
    end
