# Jevo - Evolution simulations in Julia

[![Build Status](https://travis-ci.com/tomroesch/Jevo.jl.svg?branch=master)](https://travis-ci.com/tomroesch/Jevo.jl)
[![Build Status](https://travis-ci.org/tomroesch/Jevo.jl.svg?branch=master)](https://travis-ci.org/tomroesch/Jevo.jl)

## Why?
This package enables 

## Setting up the environment

All code in this repository is written in Julia ([version 1.6](https://github.com/JuliaLang/julia/releases/tag/v1.6.0), can be installed from this link). Once Julia is installed, the standard Julia REPL can be started from the terminal, by navigating into the folder where Julia was installed, and running `path/julia/bin/julia`. On Mac, it should simply be `/Applications/Julia-1.x.app/Contents/Resources/julia/bin/julia` (we recommend setting an alias, x stands for the version of Julia that is installed).

To add this package, navigate into this folder and run: 
```julia
julia> ] dev .
```

This adds the package in the development mode, and lets you add and change things at will. Make sure to include testing for additions to the package to ensure functionality. To test the package, the environment has to the activated. Therefore, run Julia with this directory as working directory and use:

```julia
julia> ] activate .

julia> ] test
```

## Basic Use
### Initializing a Population

Included in this package are a few number of population types of genetic possible genetic sequences composed of integers as letters. Population types are initialized as subtypes of the type `populations`.

```julia
julia> using Jevo
julia>? Jevo.populations
  No documentation found.

  Summary
  ≡≡≡≡≡≡≡≡≡

  abstract type Jevo.populations <: Any

  Subtypes
  ≡≡≡≡≡≡≡≡≡≡

  Jevo.DT_population
  Jevo.simple_sites
```

The most basic type is `simple_sites`, which consists of a population of single sequences with fixed length. To find information about the fields, simply check out the documentation.

```julia
julia>?Jevo.simple_sites
  mutable struct simple_sites <: populations

  Population of simple binding sites with fixed length.

  ...

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •    N::Float64=1000: population size.

    •    l::Int64=10: Length of site.

    •    n::Int64=4: Alphabet size.

    •    seqs::Array{Array{Int64, 1}, 1}=[]: Array of sequences for each species.

    •    freqs::Array{Int64, 1}=[] Species sizes.

  ...

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> Jevo.simple_sites()
  Jevo.simple_sites(1000, 10, 4, Array{Int64,1}[], Int64[])
```

By default, populations are initiated without sequences, which are initiated in a later step. All parameters are given as keyword arguments. Initial sequences and frequencies can also be given if wanted.

```julia
julia> Jevo.simple_sites()
  Jevo.simple_sites(1000, 10, 4, Array{Int64,1}[], Int64[])

julia> Jevo.simple_sites(N=100, n=2, l=20)
Jevo.simple_sites(100, 20, 2, Array{Int64,1}[], Int64[])
```

If no initial sequences were given, they can be added by using the `initiate!` function, which takes a population and optionally the number of species as arguments (by default a single sequence is generated). Sequences are created randomly.

```julia
julia> pop = Jevo.simple_sites()
Jevo.simple_sites(1000, 10, 4, Array{Int64,1}[], Int64[])

julia> Jevo.initiate!(pop)

julia> pop
Jevo.simple_sites(1000, 10, 4, [[4, 2, 1, 4, 1, 1, 2, 2, 4, 2]], [1000])

julia> pop = Jevo.simple_sites()
Jevo.simple_sites(1000, 10, 4, Array{Int64,1}[], Int64[])

julia> Jevo.initiate!(pop, 5)

julia> pop.seqs
5-element Array{Array{Int64,1},1}:
 [2, 1, 4, 3, 2, 2, 3, 3, 2, 1]
 [2, 2, 3, 1, 4, 4, 2, 4, 4, 4]
 [3, 1, 2, 1, 1, 1, 2, 2, 4, 2]
 [2, 3, 2, 1, 4, 1, 3, 3, 4, 2]
 [2, 3, 3, 4, 1, 2, 2, 4, 4, 2]

julia> pop.freqs
5-element Array{Int64,1}:
 200
 200
 200
 200
 200
```

Other population types have more parameters, but are initialized in the same way. To see all population types that are included at the moment use the function `showtypetree`, which is copied from [wikibooks](https://en.wikibooks.org/wiki/Introducing_Julia/Types).

```
julia> showtypetree(Jevo.populations)
Jevo.populations
	Jevo.DT_population
		Jevo.driver_trailer
		Jevo.driver_trailer_l
		Jevo.mono_pop
	Jevo.simple_sites
```

Each of these types has its own `initiate!` function, but works the same for all types. What parameters are used for each type can be found in the documentation of each type. If a custom population is added to the package and used, make sure to define an appropriate `initiate!` function to the `initiate.jl` sub-module.

One can also generate a sequence which is optimally adapted, meaning it has minimal binding energy. Therefore, the second argument when initializing the sequence is an energy matrix.

