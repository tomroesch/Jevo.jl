abstract type fitness_functions end


"""
    mutable struct exponential_fitness <: fitness_functions

Fitness landscape that has the form of an exponential function.
"""
mutable struct exponential_fitness <: fitness_functions
    l::Int64
    beta::Float64
    epsilon::Float64
    f0::Float64
    fl::Float64
    l0::Int64
    E_Star::Function
end

# Setting some defaults
exponential_fitness(;l=10, beta=1, epsilon=2, f0=1, fl=0, l0=5, E_Star=Est) = exponential_fitness(l, beta, epsilon, f0, fl, l0, E_Star)


"""
    function fitness(E::Float64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Real, p::exponential_fitness)::Float64
    fitness = p.f0 * (1 - exp(p.beta * (E - p.E_Star(p.l, p.epsilon, p.l0)))) - p.fl * p.l
    return fitness
end


"""
    function fitness(E::Float64, l::Int64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Real, l::Int64, p::exponential_fitness)::Float64
    fitness = p.f0 * (1 - exp(p.beta * (E - p.E_Star(l, p.epsilon, p.l0)))) - p.fl * l
    return fitness
end



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
    epsilon::Float64
    f0::Float64
    fl::Float64
    l0::Int64
    E_Star::Function
end

# Setting some defaults
fermi_fitness(;l=10, beta=1, epsilon=2, f0=1, fl=0, l0=10, E_Star=Est) = fermi_fitness(l, beta, epsilon, f0, fl, l0, E_Star)


"""
    function fitness(E::Float64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Real, p::fermi_fitness)::Float64
    fitness = (p.f0 * (1/(1 + exp(p.beta * (E - p.E_Star(p.l, p.epsilon, p.l0))))) - p.fl * p.l)
    return fitness
end


"""
    function fitness(E::Float64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _dE_fitness(E::Real, p::fermi_fitness)
    dfitness = -p.f0 * 1/(1 + exp(p.beta * (E - p.E_Star(p.l, p.epsilon, p.l0))))^2 * p.beta * exp(p.beta * (E - p.E_Star(p.l, p.epsilon, p.l0)))
    return dfitness
end


"""
    function fitness(E::Float64, l::Int64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Real, l::Int64, p::fermi_fitness)::Float64
    fitness = (p.f0 * (1/(1 + exp(p.beta * (E - p.E_Star(l, p.epsilon, p.l0))))) - p.fl * l)
    return fitness
end


"""
    function fitness(E::Float64, p::fermi_fitness)

Evaluate fermi fitness function.
"""
function _dE_fitness(E::Real, l::Int64, p::fermi_fitness)
    dfitness = -p.f0 * 1/(1 + exp(p.beta * (E - p.E_Star(l, p.epsilon, p.l0))))^2 * p.beta * exp(p.beta * (E - p.E_Star(l, p.epsilon, p.l0)))
    return dfitness
end


"""
    Est(l)

Threshold of the fermi landscape.
"""
function Est(l::Int, epsilon::Real, l0::Int)::Real
    return (3l/4 - l0) * epsilon
end



mutable struct quadratic_fitness <: fitness_functions
    l::Int64
    c::Real
    l0::Int64
    epsilon::Real
    E_Star::Function
end

# Setting some defaults
quadratic_fitness(;l=10, c=1, l0=10, epsilon=1, E_Star=Est) = quadratic_fitness(l, c, l0, epsilon, E_Star)


"""
    function fitness(E::Float64, p::quadratic_fitness)

Evaluate quadratic fitness function.
"""
function _fitness(E::Float64, p::quadratic_fitness)
    fitness = -p.c * (E - p.E_Star(p.l, p.epsilon, p.l0))^2
    return fitness
end




"""
    function fitness(E::Float64, l::Int64, p::quadratic_fitness)

Evaluate fermi fitness function.
"""
function _fitness(E::Float64, l::Int64, p::quadratic_fitness)
    fitness = -p.c * (E - p.E_Star(l, p.epsilon, p.l0))^2
    return fitness
end


"""
    mutable struct linear_fitness <: fitness_functions

Fitness landscape that is linear.
"""
mutable struct linear_fitness <: Jevo.fitness_functions
    l::Int64
    s::Real
    fl::Real
    n::Int64
    epsilon::Real
    l0::Int64
end

# Setting some defaults
linear_fitness(;l=10, s=1, fl=0, n=4, epsilon=1, l0=10) = linear_fitness(l, s, fl, n, epsilon, l0)

    
function _fitness(E::Real, p::linear_fitness)::Float64
    fitness = -p.fl * p.l + p.s/(p.l * p.epsilon) * ((p.n-1)/p.n * (p.l * p.epsilon)  - p.l0  - E )
    return fitness
end


function _fitness(E::Real, l::Int64, p::linear_fitness)::Float64
    fitness = -p.fl * l + p.s/(p.l * p.epsilon) * ((p.n-1)/p.n * (p.l * p.epsilon)  - p.l0 - E )
    return fitness
end

#=
"""
    mutable struct linear_fitness <: fitness_functions

Fitness landscape that has the form of a fermi function.
"""
mutable struct semi_linear_fitness <: Jevo.fitness_functions
    l::Int64
    f0::Real
    fl::Float64
    n::Int64
    epsilon::Float64
end

# Setting some defaults
semi_linear_fitness(;l=10, f0=1, fl=0, n=4, epsilon=1) = semi_linear_fitness(l, f0, fl, n, epsilon)

    
function _fitness(E::Real, p::semi_linear_fitness)::Float64
    if E < p.epsilon * ((p.n-1)/p.n *p.l  - 5)
        fitness = -p.fl * p.l + p.f0
    else 
        fitness = -p.fl * p.l - p.f0 / (5p.epsilon) * (E - p.l * p.epsilon * (p.n - 1)/p.n)
    end
    return fitness
end


function _fitness(E::Real, l::Int64, p::semi_linear_fitness)::Float64
    if E < p.epsilon * ((p.n-1)/p.n * l  - 5)
        fitness = -p.fl * l + p.f0
    else 
        fitness = -p.fl * l - p.f0 / (5p.epsilon) * (E - l * p.epsilon * (p.n - 1)/p.n)
    end
    return fitness
end
=#

# Add some beautiful broadcasting
fitness(E, l, p) = _fitness(E, l, p)
Broadcast.broadcasted(::typeof(fitness), E, l, p) = broadcast(_fitness, E, l, Ref(p))


fitness(E, p) = _fitness(E, p)
Broadcast.broadcasted(::typeof(fitness), E, p) = broadcast(_fitness, E, Ref(p))


dE_fitness(E, l, p) = _dE_fitness(E, l, p)
Broadcast.broadcasted(::typeof(dE_fitness), E, l, p) = broadcast(_dE_fitness, E, l, Ref(p))


dE_fitness(E, p) = _dE_fitness(E, p)
Broadcast.broadcasted(::typeof(dE_fitness), E, p) = broadcast(_dE_fitness, E, Ref(p))


mutable struct num_fermi <: Jevo.fitness_functions
    n::Real
    l_0::Real
    gap::Real
    f0::Real
    fl::Real
end

num_fermi(;n=4, l_0=10, gap=7, f0=400, fl=10) = num_fermi(n, l_0, gap, f0, fl)

# Neutral Expectation
γ_0(n) = (n-1)/n

# Binding Threshold
γ_1(l, n, l_0) = γ_0(n) - l_0/l

# Binding probability
pb(γ, l, n, l_0, gap) = 1 / (1 + exp(gap * (l / l_0) * (γ - γ_1(l, n, l_0))))

# Fitness component for functional binding
F_b(γ, l, n, l_0, gap, f0) = f0 * pb(γ, l, n, l_0, gap)

# Fitness component of genomic constraint
F_c(l, l_0, fl) = - fl * l / l_0

# Total Fitness
F(γ, l, n, l_0, gap, f0, fl) = F_b(γ, l, n, l_0, gap, f0) + F_c(l, l_0, fl)

Jevo._fitness(γ::Real, l, p::num_fermi)::Float64 = F(γ/(l * p.gap/p.l_0), l, p.n, p.l_0, p.gap, p.f0, p.fl) 


mutable struct const_fitness <: Jevo.fitness_functions
    s::Real
end

const_fitness(;s=1) = const_fitness(s)

function _fitness(p::const_fitness)::Real
    return p.s
end

function _fitness(E::Real, p::const_fitness)::Real
    return p.s
end

function _fitness(E::Real, l::Int64, p::const_fitness)::Real
    return p.s
end

fitness(p::const_fitness) = _fitness(p::const_fitness)