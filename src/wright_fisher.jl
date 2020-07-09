using Distributions

function get_energy(pop::binding_sites, emat::Array{T, 2}) where {T<:Real}
    c = length(pop.seqs)
    energy = zeros(Float64, c)
    for j in 1:c
        for i in eachindex(1:pop.l)
            energy[j] += emat[pop.seqs[j][i], i]
        end
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
    function sample!(pop::binding_sites, fit::fitness_functions, emat::Array{T, 2}) where {T<:Real}

Sample a new generation and remove extinct species.
"""
function sample_gen!(pop::binding_sites, fit::fitness_functions, emat::Array{T, 2}) where {T<:Real}
    E = get_energy(pop, emat)
    f::Array{Float64, 1} = fitness.(E, fit)
    mean_fitness::Float64 = sum(1 / pop.N * pop.freqs.* f)
    norm::Float64 = sum(pop.freqs .* exp.(f .- mean_fitness))
    probabilities::Array{Float64, 1} = pop.freqs .* exp.(f .- mean_fitness) ./ norm
    pop.freqs = rand(Multinomial(pop.N, probabilities), 1)[:]
    remove_empty!(pop.freqs, pop.seqs)
end