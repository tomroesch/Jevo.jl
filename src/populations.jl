abstract type populations end


mutable struct binding_sites <: populations
    # Population size
    N::Int64
    # Maximal length (in case of length mutations)
    L::Int64
    # Inital length of specific site
    l::Int64
    # Alphabet size
    n::Int64
    # Sequences
    seqs::Array{Array{Int64, 1}, 1}
    # Subpopulation sizes
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


mutable struct driver_trailer <: populations
    # Population size
    N::Int64
    # Maximal length (in case of length mutations)
    L::Int64
    # Inital length of specific site
    l::Int64
    # Driver alphabet size
    m::Int64
    # Trailer alphabet size
    n::Int64
    # Driver sequence
    driver::Array{Int64, 1}
    # Trailer sequences
    seqs::Array{Array{Int64, 1}, 1}
    # Subpopulation sizes
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