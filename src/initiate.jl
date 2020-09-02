export initiate!,make_driver

"""
    function initiate!(pop::binding_sites, c::Int64; overwrite=false)

Iniate a population of random sequences.

# Arguments
- `pop::binding_sites`: population that is to be filled.
- `c::Int64=1`: number of species to create.
- `overwrite::Bool=false`: If `true`, resets a given population to equal subspecies size.
"""
function initiate!(pop::binding_sites, c::Int64=1; overwrite::Bool=false)

    if ~isempty(pop.seqs)
        # Reiniate existing sequences
        if overwrite
            c = length(pop.seqs)
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
    function initiate!(
        pop::driver_trailer,
        c::Int64;
        driver::Array{Int64, 1}=Int64[],
        overwrite::Bool=false,
        opt::Bool=false)

Initiates population.


Creates a population of type `driver_trailer`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.

# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `c::Int64=1`: number of species to create.
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
- `overwrite::Bool=false`: If `true`, resets a given population to equal subspecies size.
- `opt::Bool="false"`: If `true`, copies driver to seq. If `false`, creates random sequence.
"""
function initiate!(
    pop::driver_trailer,
    c::Int64=1;
    driver::Array{Int64, 1}=Int64[],
    overwrite::Bool=false,
    opt::Bool=false
    )

    if ~isempty(pop.seqs)
        # Reiniate existing sequences
        if overwrite
            c = length(pop.seqs)
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
        pop.driver[:] = make_driver(driver, pop.m, pop.L)
        # Create new population
        if ~opt
            N_sub = pop.N ÷ c
            rest = pop.N - N_sub * c
            for i in 1:c
                push!(pop.seqs, rand(collect(1:pop.n), pop.L))
                push!(pop.freqs, N_sub)
            end
            pop.freqs[1] += rest
        else
            push!(pop.freqs, pop.N)
            push!(pop.seqs, pop.driver)
        end
    end
end


"""
    function initiate!(
        pop::driver_trailer_l,
        c::Int64;
        driver::Array{Int64, 1}=Int64[],
        overwrite=false,
        opt::Bool=false)

Iniate a population of c species of and a driver sequence, which can be given.

Creates a population of type `driver_trailerPl`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.


# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `c::Int64=1`: number of species to create.
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
- `overwrite::Bool=false`: If `true`, resets a given population to equal subspecies size.
- `opt::Bool="false"`: If `true`, copies driver to seq. If `false`, creates `c` random sequences.
"""
function initiate!(
    pop::driver_trailer_l,
    c::Int64=1;
    driver::Array{Int64, 1}=Int64[],
    overwrite::Bool=false,
    opt::Bool=false
    )

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
            pop.l = ones(Int64, c) * pop.l_0
        else
            return
        end
    else
        pop.driver[:] = make_driver(driver, pop.m, pop.L)

        if ~opt
            # Create new population
            N_sub = pop.N ÷ c
            rest = pop.N - N_sub * c
            for i in 1:c
                push!(pop.seqs, rand(collect(1:pop.n), pop.L))
                push!(pop.freqs, N_sub)
            end
            pop.freqs[1] += rest

            pop.l = ones(Int64, c) * pop.l_0
        else
            push!(pop.seqs, pop.driver)
            push!(pop.freqs, pop.N)
            push!(pop.l, pop.l_0)
        end

    end
end



"""
    function initiate!(
        pop::mono_pop;
        driver::Array{Int64, 1}=Int64[],
        opt::Bool=false)

Iniate a monomorphic population and a driver sequence, which can be given.

Creates a population of type `mono_pop`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.


# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
- `opt::Bool=false`: If `false`, creates random sequence. If `true`, copies driver sequence.
"""
function initiate!(
    pop::mono_pop;
    driver::Array{Int64, 1}=Int64[],
    opt::Bool=false
    )

    pop.driver = make_driver(driver, pop.m, pop.l)

    if ~opt
        pop.seqs =  rand(collect(1:pop.n), pop.l)
    else
        pop.seqs = pop.driver
    end
end


"""
    make_driver(given_driver, m, L)

Create a driver sequence.
"""
function make_driver(given_driver, m, L)
    driver = zeros(Int64, L)
    if isempty(given_driver)
        driver[:] = rand(collect(1:m), L)
    else
        if any(given_driver .> m)
            throw(ArgumentError("Driver sequence has elements outside of alphabet."))
        elseif any(given_driver .< 0)
            throw(ArgumentError("Don't include negative elements in driver sequence."))
        end
        if length(given_driver) == L
            driver[:] = given_driver
        elseif length(given_driver) < L
            driver[1:length(given_driver)] = given_driver
            driver[length(given_driver) + 1:end] = rand(collect(1:m), L - length(given_driver))
        else
            driver[:] = given_driver[1:L]
        end
    end
    return driver
end
