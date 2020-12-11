# Jevo - Evolution simulations in Julia

[![Build Status](https://travis-ci.com/tomroesch/Jevo.jl.svg?branch=master)](https://travis-ci.com/tomroesch/Jevo.jl)
[![Build Status](https://travis-ci.org/tomroesch/Jevo.jl.svg?branch=master)](https://travis-ci.org/tomroesch/Jevo.jl)

### Setting up the environment

All code in this repository is written in Julia ([version 1.5](https://github.com/JuliaLang/julia/releases/tag/v1.5.0), can be installed from this link). Once Julia is installed, the standard Julia REPL can be started from the terminal, by navigating into the folder where Julia was installed, and running `path/julia/bin/julia`. On Mac, it should simply be `/Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia` (we recommend setting an alias).

To add this package, navigate into this folder and run 
```
julia> ] dev .
```
which adds the package in the development mode, and lets you add and change things at will. Make sure to include testing for additions to the package to ensure functionality.