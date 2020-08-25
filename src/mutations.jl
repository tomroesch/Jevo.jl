using StatsBase

export mutation!,driver_mutation!,length_mutation!

"""
    function mutation!(pop::populations)

Mutate a base pair in a random species, chosen depending on subpopulation size.
If the species is extinct after mutation, it is removed. If mutation creates a new species,
it is added to the population, otherwise the subpopulation size of species the individual mutates to
is simply increased by one.
"""
function mutation!(pop::populations)
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
    function mutation!(pop::driver_trailer_l)

Mutate a base pair in a random species, chosen depending on subpopulation size.
If the species is extinct after mutation, it is removed. If mutation creates a new species,
it is added to the population, otherwise the subpopulation size of species the individual mutates to
is simply increased by one.
"""
function mutation!(pop::driver_trailer_l)
    # Choose random sequence depending on subpopulation size
    mutation_sequence = sample(collect(1:length(pop.seqs)), Weights(pop.freqs ./ pop.N))
    # Add sequence to mutate
    push!(pop.seqs, deepcopy(pop.seqs[mutation_sequence]))
    push!(pop.freqs, 1)
    push!(pop.l, pop.l[mutation_sequence])
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
        deleteat!(pop.l, length(pop.l))
    end
    # Remove extinct species
    remove_empty!(pop.freqs, pop.seqs, pop.l)
end

"""
    function driver_mutation!(pop::DT_population)

Mutate a random base in the driver sequence.
"""
function driver_mutation!(pop::DT_population)
    # Choose random base to mutate
    mutation_base = rand(1:pop.L)
    # Exclude self mutations
    poss_bases = filter(x->x != pop.driver[mutation_base], collect(1:pop.m))
    # Choose random new base
    pop.driver[mutation_base] = rand(poss_bases)
end


"""
    function driver_mutation!(pop::driver_trailer)

Mutate the length of species and create a new species
"""
function length_mutation!(pop::driver_trailer_l)
    # Choose random sequence depending on subpopulation size
    mutation_sequence = sample( collect(1:length(pop.seqs)),
                                Weights(pop.freqs ./ pop.N))
    # Add sequence to mutate

    if rand() < 0.5 && pop.l[mutation_sequence] < pop.L
        var = 1
    elseif pop.l[mutation_sequence] > 2
        var = -1
    else
        return
    end

    push!(pop.seqs, copy(pop.seqs[mutation_sequence]))
    push!(pop.freqs, 1)
    push!(pop.l, pop.l[mutation_sequence])

    pop.freqs[mutation_sequence] -= 1
    pop.l[end] += var
    # Remove extinct species
    remove_empty!(pop.freqs, pop.seqs, pop.l)
end


"""
    bp_substitution!(pop::populations, emat::Array{T, 2}, fitness_function::fitness_functions) where {T<:Real}

Attempt a base pair substitution. The population has to consist of one species
only. Mutation is accepted with probability given by the Kimura fixation
probability.
"""
function bp_substitution!(pop::populations, emat::Array{T, 2}, fitness_function::fitness_functions) where {T<:Real}
    if length(pop.seqs) > 1
        throw(ArgumentError("Input population has more than one species."))
    end
    # Choose random base to mutate
    mutation_base = rand(1:pop.L)
    # Exclude self mutations
    poss_bases = filter(x->x != pop.seqs[1][mutation_base], collect(1:pop.n))
    # Choose random new base
    temp_pop = deepcopy(pop)
    temp_pop.seqs[1][mutation_base] = rand(poss_bases)
    # Compute energies
    E = get_energy(pop, emat)
    E_mutant = get_energy(temp_pop, emat)
    # Compute fitness and selection coefficients
    F = fitness(E[1], pop.l[1], fitness_function)
    F_mutant = fitness(E_mutant[1], pop.l[1], fitness_function)
    s = F_mutant - F
    # Accept mutatant depending on Kimura probability
    if rand() < kimura_prob(s, pop.N)
        pop.seqs = temp_pop.seqs
    end
end


"""
    l_substitution!(
        pop::driver_trailer_l,
        emat::Array{T, 2},
        fitness_function::fitness_functions) where {T<:Real}

Attempt a substitution of a length mutation. Population has to be consisting
of one species only.
"""
function l_substitution!(
    pop::driver_trailer_l,
    emat::Array{T, 2},
    fitness_function::fitness_functions) where {T<:Real}

    if length(pop.seqs) > 1
        throw(ArgumentError("Input population has more than one species."))
    end

    if rand() < 0.5 && pop.l[1] < pop.L
        var = 1
    elseif pop.l[1] > 2
        var = -1
    else
        return
    end
    # Choose random new base
    temp_pop = deepcopy(pop)
    temp_pop.l[1] += var
    # Compute energies
    E = get_energy(pop, emat)
    E_mutant = get_energy(temp_pop, emat)
    # Compute fitness and selection coefficients
    F = fitness(E[1], pop.l[1], fitness_function)
    F_mutant = fitness(E_mutant[1], temp_pop.l[1], fitness_function)
    s = F_mutant - F
    # Accept mutatant depending on Kimura probability
    if rand() < kimura_prob(s, pop.N)
        pop.l = temp_pop.l
    end
end


"""
    kimura_prob(s, N)

Probability of fixation of a mutation in a monomorphic population.
"""
function kimura_prob(s, N)
    if abs(s) <= 10^-8
        return 1 / N
    else
        return (1 - exp(-2s)) / (1 - exp(-2N * s))
    end
end
