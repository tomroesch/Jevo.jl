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
    freqs::Array{Int64}
end

binding_sites(;N=10000, L=10, l=10, n=4, seqs=Array{Int64, 1}[], freqs=Int64[]) = binding_sites(N, L, l, n, seqs, freqs)

    