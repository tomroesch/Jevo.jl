export fermi_fitness, fitness

abstract type fitness_functions end


"""
    mutable struct fermi_fitness <: fitness_functions

Fitness landscape that has the form of a fermi function.

# Defaults

```
fermi_fitness(;l=l, beta=beta, f0=f0) = fermi_fitness(l, beta, f0, 0, Est)
fermi_fitness(;beta=beta, f0=f0) = fermi_fitness(10, beta, f0, 0, Est)
fermi_fitness(;l=l, beta=beta, f0=f0, fl=fl) = fermi_fitness(l, beta, f0, fl, Est)
```
"""
mutable struct fermi_fitness <: fitness_functions
    l::Int64
    beta::Float64
    f0::Float64
    fl::Float64
    E_Star::Function
end

# Setting some defaults
fermi_fitness(;l=l, beta=beta, f0=f0) = fermi_fitness(l, beta, f0, 0, Est)
fermi_fitness(;beta=beta, f0=f0) = fermi_fitness(10, beta, f0, 0, Est)
fermi_fitness(;l=l, beta=beta, f0=f0, fl=fl) = fermi_fitness(l, beta, f0, fl, Est)


"""
    function fitness(E::Float64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Float64, p::fermi_fitness)
    fitness = (p.f0 * (1/(1 + exp(-p.beta * (E - p.E_Star(p.l))))) - p.fl * p.l)
    return fitness
end

# Add some beautiful broadcasting
fitness(E, p) = _fitness(E, p)
Broadcast.broadcasted(::typeof(fitness), E, p) = broadcast(_fitness, E, Ref(p))


"""
    function fitness(E::Float64, l::Int64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Float64, l::Int64, p::fermi_fitness)
    fitness = (p.f0 * (1/(1 + exp(-p.beta * (E - p.E_Star(l))))) - p.fl * l)
    return fitness
end


# Add some beautiful broadcasting
fitness(E, l, p) = _fitness(E, l, p)
Broadcast.broadcasted(::typeof(fitness), E, l, p) = broadcast(_fitness, E, l, Ref(p))

"""
    Est(l)

Threshold of the fermi landscape.
"""
function Est(l)
    return (l/4 + 4)
end


