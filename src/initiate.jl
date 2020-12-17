
"""
    function initiate!(pop::binding_sites, c::Int64; overwrite=false)

Iniate a population of random sequences.

# Arguments
- `pop::binding_sites`: population that is to be filled.
- `c::Int64=1`: number of species to create.
"""
function initiate!(pop::binding_sites, c::Int64=1)
    # Create new population
    N_sub = pop.N ÷ c
    rest = pop.N - N_sub * c
    for i in 1:c
        push!(pop.seqs, rand(collect(1:pop.n), pop.l))
        push!(pop.freqs, N_sub)
    end
    pop.freqs[1] += rest
    nothing
end



"""
    function initiate!(
        pop::driver_trailer,
        c::Int64;
        driver::Array{Int64, 1}=Int64[]
        )

Initiates population.


Creates a population of type `driver_trailer`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.

# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `c::Int64=1`: number of species to create.
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
"""
function initiate!(
    pop::driver_trailer,
    c::Int64=1;
    driver::Array{Int64, 1}=Int64[]
    )

    pop.driver[:] = make_driver(driver, pop.m, pop.L)
    # Create new population
    N_sub = pop.N ÷ c
    rest = pop.N - N_sub * c
    for i in 1:c
        push!(pop.seqs, rand(collect(1:pop.n), pop.L))
        push!(pop.freqs, N_sub)
    end
    pop.freqs[1] += rest
    nothing
end


"""
    function initiate!(
        pop::driver_trailer,
        emat::Array{Float64, 2};
        driver::Array{Int64, 1}=Int64[]
        )

Initiates an optimal population.


Creates a population of type `driver_trailer`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Optimal sequence
will chosen as the one which minimizes energy.

# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `emat::Array{Float64,2}` : Energy matrix to choose optimal sequence
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
"""
function initiate!(
    pop::driver_trailer,
    emat::Array{Float64, 2};
    driver::Array{Int64, 1}=Int64[]
    )

    if size(emat) != (pop.n, pop.m)
        throw(DimensionMismatch("Energy matrix needs to have dimensions (n, m)"))
    end

    pop.driver[:] = make_driver(driver, pop.m, pop.L)
    push!(pop.freqs, pop.N)

    push!(pop.seqs, argmin.([emat[:, j] for j in pop.driver]))

end



"""
    function initiate!(
        pop::driver_trailer_l,
        c::Int64;
        driver::Array{Int64, 1}=Int64[]
        )

Iniate a population of c species of and a driver sequence, which can be given.

Creates a population of type `driver_trailerPl`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.


# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `c::Int64=1`: number of species to create.
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
"""
function initiate!(
    pop::driver_trailer_l,
    c::Int64=1;
    driver::Array{Int64, 1}=Int64[],
    )

    pop.driver[:] = make_driver(driver, pop.m, pop.L)
    N_sub = pop.N ÷ c
    rest = pop.N - N_sub * c
    for i in 1:c
        push!(pop.seqs, rand(collect(1:pop.n), pop.L))
        push!(pop.freqs, N_sub)
    end
    pop.freqs[1] += rest

    pop.l = ones(Int64, c) * pop.l_0
end


"""
    function initiate!(
        pop::driver_trailer_l,
        emat::Array{Float64, 2};
        driver::Array{Int64, 1}=Int64[]
        )

Iniate a population of c species of and a driver sequence, which can be given.

Creates a population of type `driver_trailerPl`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.


# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `emat::Array{Float64,2}` : Energy matrix to choose optimal sequence
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
"""
function initiate!(
    pop::driver_trailer_l,
    emat::Array{Float64,2};
    driver::Array{Int64, 1}=Int64[],
    )

    pop.driver[:] = make_driver(driver, pop.m, pop.L)
    push!(pop.freqs, pop.N)
    push!(pop.seqs, argmin.([emat[:, j] for j in pop.driver]))
        
    pop.l = Int64[pop.l_0]
end


"""
    function initiate!(
        pop::mono_pop;
        driver::Array{Int64, 1}=Int64[]
        )

Iniate a monomorphic population and a driver sequence, which can be given.

Creates a population of type `mono_pop`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.


# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
"""
function initiate!(
    pop::mono_pop;
    driver::Array{Int64, 1}=Int64[],
    )
    pop.driver = make_driver(driver, pop.m, pop.l)
    pop.seqs =  rand(collect(1:pop.n), pop.l)
    nothing
end


"""
    function initiate!(
        pop::mono_pop
        emat::Array{Float64, 2};
        driver::Array{Int64, 1}=Int64[]
        )

Iniate a monomorphic population and a driver sequence, which can be given.

Creates a population of type `mono_pop`, where a driver sequence can be given.
Given driver sequence can be shorter or longer than initial binding site length. Sequence
will be either filled with random letters or truncated to the right length.


# Arguments
- `pop::driver_trailer`: population that is to be filled.
- `emat::Array{Float64,2}` : Energy matrix to choose optimal sequence
- `driver::Array{Int64, 1}=Int64[]`: Optional initial driver sequence.
"""
function initiate!(
    pop::mono_pop,
    emat::Array{Float64,2};
    driver::Array{Int64, 1}=Int64[],
    )
    if size(emat) != (pop.n, pop.m)
        throw(DimensionMismatch("Energy matrix needs to have dimensions (n, m)"))
    end
    pop.driver = make_driver(driver, pop.m, pop.l)
    pop.seqs =  argmin.([emat[:, j] for j in pop.driver])
    nothing
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
