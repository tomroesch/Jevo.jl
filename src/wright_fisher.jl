using Distributions


"""
    get_energy(pop::binding_sites, emat::Array{T, 2}) where {T<:Real}

Compute binding energies in simple binding energy model. Energy matrix needs to have
dimensions (n, site_length).
"""
function get_energy(pop::binding_sites, emat::Array{T, 2}) where {T<:Real}
    c = length(pop.seqs)
    energy = zeros(T, c)
    for j in 1:c
        for i in eachindex(1:pop.l)
            energy[j] += emat[pop.seqs[j][i], i]
        end
    end
    return energy
end


"""
    get_energy(pop::driver_trailer, emat::Array{T, 2}) where {T<:Real}

Compute binding energies in driver/trailer population.
"""
function get_energy(pop::driver_trailer, emat::Array{T, 2}) where {T<:Real}
    c = length(pop.seqs)
    energy = zeros(Float64, c)
    for j in 1:c
        for i in eachindex(1:pop.l)
            energy[j] += emat[pop.seqs[j][i], pop.driver[i]]
        end
    end
    return energy
end


"""
    get_energy(pop::driver_trailer_l, emat::Array{T, 2}) where {T<:Real}

Compute binding energies in driver/trailer population with different binding lengths.
"""
function get_energy(pop::driver_trailer_l, emat::Array{T, 2}) where {T<:Real}
    c = length(pop.seqs)
    energy = zeros(Float64, c)
    for j in 1:c
        for i in 1:pop.l[j]
            energy[j] += emat[pop.seqs[j][i], pop.driver[i]]
        end
    end
    return energy
end


"""
    get_energy(pop::mono_pop, emat::Array{T, 2}) where {T<:Real}

Compute binding energies in monomorphic population..
"""
function get_energy(pop::mono_pop, emat::Array{T, 2})::T where {T<:Real}

    energy::T = 0
    for i in eachindex(pop.seqs)
        energy += emat[pop.seqs[i], pop.driver[i]]
    end
    return energy
end


function remove_empty!(array_list...)
    main_arr = array_list[1]
    inds = findall(x->x==0, main_arr)
    for array in array_list
        deleteat!(array, inds)
    end
end


"""
    function sample_gen!(pop::binding_sites, fit::fitness_functions, emat::Array{T, 2}; remove=true) where {T<:Real}

Sample a new generation and remove extinct species.
"""
function sample_gen!(pop::populations, fit::fitness_functions, emat::Array{T, 2}; remove=true) where {T<:Real}
    E = get_energy(pop, emat)
    f::Array{Float64, 1} = fitness.(E, fit)
    mean_fitness::Float64 = sum(1 / pop.N * pop.freqs.* f)
    norm::Float64 = sum(pop.freqs .* exp.(f .- mean_fitness))
    probabilities::Array{Float64, 1} = pop.freqs .* exp.(f .- mean_fitness) ./ norm
    pop.freqs = rand(Multinomial(pop.N, probabilities), 1)[:]
    if remove
        remove_empty!(pop.freqs, pop.seqs)
    end
end



"""
    function sample_gen!(pop::driver_trailer, fit::fitness_functions, emat::Array{T, 2}; remove=true) where {T<:Real}

Sample a new generation and remove extinct species.
"""
function sample_gen!(pop::driver_trailer_l, fit::fitness_functions, emat::Array{T, 2}; remove=true) where {T<:Real}
    E = get_energy(pop, emat)
    f::Array{Float64, 1} = fitness.(E, pop.l, fit)
    mean_fitness::Float64 = sum(1 / pop.N * pop.freqs.* f)
    norm::Float64 = sum(pop.freqs .* exp.(f .- mean_fitness))
    probabilities::Array{Float64, 1} = pop.freqs .* exp.(f .- mean_fitness) ./ norm
    pop.freqs = rand(Multinomial(pop.N, probabilities), 1)[:]
    if remove
        remove_empty!(pop.freqs, pop.seqs, pop.l)
    end
end
