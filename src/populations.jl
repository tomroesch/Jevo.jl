abstract type populations end

abstract type DT_population <: populations end

"""
    mutable struct binding_sites <: populations

Creates population of simple binding sites of fixed length.
"""
mutable struct binding_sites <: populations
    "Population size"
    N::Int64
    "Maximal length (in case of length mutations)"
    L::Int64
    "Inital length of specific site"
    l::Int64
    "Alphabet size"
    n::Int64
    "Sequences"
    seqs::Array{Array{Int64, 1}, 1}
    "Subpopulation sizes"
    freqs::Array{Int64, 1}

    function binding_sites(N, L, l, n, seqs, freqs)
        if l > L
            @warn "l>L. Choosing L=l."
            L = l
        end
        new(N, L, l, n, seqs, freqs)
    end
end

binding_sites(;N=10000, L=10, l=L, n=4, seqs=Array{Int64, 1}[], freqs=Int64[]) = binding_sites(N, L, l, n, seqs, freqs)


"""
    mutable struct driver_trailer <: DT_population

Creates population of simple binding sites of fixed length with a given driver sequence.
"""
mutable struct driver_trailer <: DT_population
    "Population size"
    N::Int64
    "Maximal length (in case of length mutations)"
    L::Int64
    "Inital length of specific site"
    l::Int64
    "Driver alphabet size"
    m::Int64
    "Trailer alphabet size"
    n::Int64
    "Driver sequence"
    driver::Array{Int64, 1}
    "Trailer sequences"
    seqs::Array{Array{Int64, 1}, 1}
    "Subpopulation sizes"
    freqs::Array{Int64, 1}

    function driver_trailer(N, L, l, m, n, driver, seqs, freqs)
        if l > L
            @warn "l>L. Choosing L=l."
            L = l
        end
        new(N, L, l, m, n, driver, seqs, freqs)
    end
end

driver_trailer(;N=10000, L=10, l=L, m=4, n=4, driver=zeros(Int64, L), seqs=Array{Int64, 1}[], freqs=Int64[]) = driver_trailer(N, L, l, m, n, driver, seqs, freqs)


"""
    mutable struct driver_trailer_l <: DT_population

Creates population of binding sites with a common driver sequence, variable length.
"""
mutable struct driver_trailer_l <: DT_population
    "Population size"
    N::Int64
    "Maximal length (in case of length mutations)"
    L::Int64
    "Inital length of specific site"
    l_0::Int64
    "Driver alphabet size"
    m::Int64
    "Trailer alphabet size"
    n::Int64
    "Driver sequence"
    driver::Array{Int64, 1}
    "Trailer sequences"
    seqs::Array{Array{Int64, 1}, 1}
    "Subpopulation sizes"
    freqs::Array{Int64, 1}
    "length for each subpopulation"
    l::Array{Int64, 1}

    function driver_trailer_l(N, L, l_0, m, n, driver, seqs, freqs, l)
        if l_0 > L
            @warn "l_0>L. Choosing L=l_0."
            L = l_0
        end
        new(N, L, l_0, m, n, driver, seqs, freqs, l)
    end
end

driver_trailer_l(;N=10000, L=10, l_0=L, m=4, n=4, driver=zeros(Int64, L), seqs=Array{Int64, 1}[], freqs=Int64[], l=Int64[]) = driver_trailer_l(N, L, l_0, m, n, driver, seqs, freqs, l)


"""
    mutable struct mono_pop <: populations

Creates population of a single binding site and driver.
"""
mutable struct mono_pop <: DT_population
    "Population size"
    N::Int64
    "Inital length of specific site"
    l::Int64
    "Alphabet size"
    n::Int64
    "Sequences"
    seqs::Array{Int64, 1}
    "Driver alphabet size"
    m::Int64
    "Driver sequence"
    driver::Array{Int64, 1}
end

mono_pop(;N=1000, l=15, n=4, seqs=Int64[], m=4, driver=Int64[]) = mono_pop(N, l, n, seqs, m, driver)
