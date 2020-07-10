"""
    function initiate_rand!(pop::binding_sites, c::Int64; overwrite=false)

Iniate a population of c species of random sequences.
"""
function initiate_rand!(pop::binding_sites, c::Int64; overwrite=false)
    
    if ~isempty(pop.seqs)
        # Reiniate existing sequences 
        if overwrite
            c = size(pops.seqs)
            pop.freqs .= 0
            N_sub = pop.N ÷ c
            rest = pop.N - N_sub * c
            for i in 1:c
                pop.freqs[i] =  N_sub
            end
            pop.freqs[1] += rest
        else
            return
        end
    else
        # Create new population
        N_sub = pop.N ÷ c
        rest = pop.N - N_sub * c
        for i in 1:c
            push!(pop.seqs, rand(collect(1:pop.n), pop.L))
            push!(pop.freqs, N_sub)
        end
        pop.freqs[1] += rest
    end
end
    


"""
    function initiate_rand!(pop::driver_trailer, c::Int64; driver::Array{Int64, 1}=Int64[], overwrite=false)

Iniate a population of c species of random sequences and a driver sequence, which can be given.
If the given sequence is shorter than the required length, the remaining positions will be filled randomly.
if the given sequence is longer than the required length, only the first `L` bases will be taken.
"""
function initiate_rand!(pop::driver_trailer, c::Int64; driver::Array{Int64, 1}=Int64[], overwrite=false)
    
    if ~isempty(pop.seqs)
        # Reiniate existing sequences 
        if overwrite
            c = size(pops.seqs)
            pop.freqs .= 0
            N_sub = pop.N ÷ c
            rest = pop.N - N_sub * c
            for i in 1:c
                pop.freqs[i] =  N_sub
            end
            pop.freqs[1] += rest
        else
            return
        end
    else
        # Create new population
        N_sub = pop.N ÷ c
        rest = pop.N - N_sub * c
        for i in 1:c
            push!(pop.seqs, rand(collect(1:pop.n), pop.L))
            push!(pop.freqs, N_sub)
        end
        pop.freqs[1] += rest
        if isempty(driver)
            pop.driver[:] = rand(collect(1:pop.m), pop.L)
        else
            if length(driver) == pop.L
                pop.driver[:] = driver
            elseif length(driver) < pop.L
                pop.driver[1:length(driver)] = driver
                pop.driver[length(driver)+1:end] = rand(collect(1:pop.m), pop.L-length(driver))
            else
                pop.driver[:] = driver[1:pop.L]
            end
        end
    end
end