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

binding_sites(;N=N, L=L, l=l, n=n) = binding_sites(N, L, l, n, Array{Int64, 1}[], Int64[])
binding_sites(;N=N, L=L, l=l) = binding_sites(N, L, l, 4, Array{Int64, 1}[], Int64[])
binding_sites(;N=N, l=l, n=n) = binding_sites(N, l, l, n, Array{Int64, 1}[], Int64[])
binding_sites(;N=N, l=l) = binding_sites(N, l, l, 4, Array{Int64, 1}[], Int64[])
    