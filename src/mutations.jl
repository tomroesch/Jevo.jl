using StatsBase


"""
    function mutation!(pop::binding_sites)

Mutate a base pair in a random species, chosen depending on subpopulation size.
If the species is extinct after mutation, it is removed. If mutation creates a new species,
it is added to the population, otherwise the subpopulation size of species the individual mutates to
is simply increased by one.
"""
function mutation!(pop::binding_sites)
    # Choose random sequence depending on subpopulation size
    mutation_sequence = sample(collect(1:length(pop.seqs)), Weights(pop.freqs ./ pop.N))
    # Add sequence to mutate
    push!(pop.seqs, copy(pop.seqs[mutation_sequence]))
    push!(pop.freqs, 1)
    pop.freqs[mutation_sequence] -= 1
    # Choose random base to mutate
    mutation_base = rand(1:pop.L)
    # Exclude self mutations
    poss_bases = filter(x->x != pop.seqs[end][mutation_base], collect(1:pop.n))
    # Choose random new base
    pop.seqs[end][mutation_base] = rand(poss_bases)
    copied = findfirst(x->~any(.~(x.==pop.seqs[end])), pop.seqs[1:end-1])
    # Look if mutated sequence already exists
    if ~isnothing(copied)
        pop.freqs[copied] += 1
        deleteat!(pop.freqs, length(pop.seqs))
        deleteat!(pop.seqs, length(pop.seqs))
    end
    # Remove extinct species
    remove_empty!(pop.freqs, pop.seqs)
end


"""
    function mutation!(pop::driver_trailer)

Mutate a base pair in a random species, chosen depending on subpopulation size.
If the species is extinct after mutation, it is removed. If mutation creates a new species,
it is added to the population, otherwise the subpopulation size of species the individual mutates to
is simply increased by one.
"""
function mutation!(pop::driver_trailer)
    # Choose random sequence depending on subpopulation size
    mutation_sequence = sample(collect(1:length(pop.seqs)), Weights(pop.freqs ./ pop.N))
    # Add sequence to mutate
    push!(pop.seqs, copy(pop.seqs[mutation_sequence]))
    push!(pop.freqs, 1)
    pop.freqs[mutation_sequence] -= 1
    # Choose random base to mutate
    mutation_base = rand(1:pop.L)
    # Exclude self mutations
    poss_bases = filter(x->x != pop.seqs[end][mutation_base], collect(1:pop.n))
    # Choose random new base
    pop.seqs[end][mutation_base] = rand(poss_bases)
    copied = findfirst(x->~any(.~(x.==pop.seqs[end])), pop.seqs[1:end-1])
    # Look if mutated sequence already exists
    if ~isnothing(copied)
        pop.freqs[copied] += 1
        deleteat!(pop.freqs, length(pop.seqs))
        deleteat!(pop.seqs, length(pop.seqs))
    end
    # Remove extinct species
    remove_empty!(pop.freqs, pop.seqs)
end


"""
    function driver_mutation!(pop::driver_trailer)

Mutate a random base in the driver sequence.
"""
function driver_mutation!(pop::driver_trailer)
    # Choose random base to mutate
    mutation_base = rand(1:pop.L)
    # Exclude self mutations
    poss_bases = filter(x->x != pop.driver[mutation_base], collect(1:pop.m))
    # Choose random new base
    pop.driver[mutation_base] = rand(poss_bases)
end
