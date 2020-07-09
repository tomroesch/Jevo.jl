"""
    createPop(P::parameters)

Create a random binding site with random preferred letters.
"""
function createPop(P::parameters)
    driver = rand(1:P.k, P.l)
    trailer = rand(1:P.n, P.l)
    return trailer,copy(trailer)
    #return trailer,driver
end


"""
    createPop(kappa::AbstractFloat, P::parameters)

Create a binding site with ratio kappa of matches.
"""
function createPop(kappa::AbstractFloat, P::parameters)
    p = pl(kappa,P)[1]
    r = rand()
    l = 1
    q = p[l]
    while q < r
        l += 1
        q += p[l]
    end
    driver = rand(1:P.k, l)
    trailer = rand(1:P.n, l)
    return trailer,driver
end


"""
    geten(trailer::Array{Int,1},driver::Array{Int,1})

Get binding energy for Berg, von Hippel model.
"""
function geten(trailer::Array{Int,1}, driver::Array{Int,1})
    energy::Float64 = 0.
    for i in eachindex(trailer)
        if driver[i] == trailer[i]
            energy += 1
        end
    end
    return energy
end


"""
    geten(trailer::Array{Int,1},driver::Array{Int,1})

Get binding energy given energy matrix.
"""
function geten(trailer::Array{Int,1},driver::Array{Int,1},enMat::AbstractArray)
    energy = sum([enMat[driver[i], trailer[i]] for i in eachindex(trailer)])
    return energy
end


"""
    tmut(trailer::Array{Int,1}, driver::Array{Int,1}, P::parameters)

Mutation in the trailer. Acception probability is given by the Kimura fraction.
"""
function tmut(trailer::Array{Int,1}, driver::Array{Int,1}, P::parameters)
    l = length(trailer)
    loci = rand(1:l)
    mutation = rand(filter(x -> x != trailer[loci], 1:P.n))
    dummie = copy(trailer)
    dummie[loci] = mutation
    s = gets(trailer, dummie, driver, P)
    prob = kimura(s, P.N)
    if prob >= rand()
        trailer[loci] = mutation
        return trailer
    else
        return trailer
    end
end


"""
    gets(trailer1::Array{Int,1}, trailer2::Array{Int,1},
            driver::Array{Int,1}, P::parameters)

Compute the selection coefficient between two genotypes.
"""
function gets(trailer1::Array{Int,1}, trailer2::Array{Int,1},
            driver::Array{Int,1}, P::parameters)

    OldE = geten(trailer1, driver)
    NewE = geten(trailer2, driver)
    OldFit = P.f(OldE, length(trailer1), P)
    NewFit = P.f(NewE, length(trailer2), P)
    s = NewFit - OldFit
    return s
end


"""
    dmut(driver::Array{Int,1}, P::parameters)

Mutation in driver sequence. Always accepted.
"""
function dmut(driver::Array{Int,1}, P::parameters)
    l = length(driver)
    amino_acid = rand(1:l)
    mutation = rand(filter(x->x!=driver[amino_acid],1:P.k))
    driver[amino_acid] = mutation
    return driver
end


function kimura(s::Float64, N)
    if abs(s) > 0.00001 && s > -2
        return  (1-exp(-2 * s)) / (1-exp(-2 * N * s))
    elseif s < -2
        return 0.
    else
        return 1 / N
    end
end


"""
    lmut(trailer::Array{Int,1}, driver::Array{Int,1}, P::parameters)

Mutation of binding site length. Acception probability given by Kimura fixation probability.
"""
function lmut(trailer::Array{Int,1}, driver::Array{Int,1}, P::parameters)
    dtrailer = deepcopy(trailer)
    ddriver = deepcopy(driver)
    if 0.5 > rand()
        push!(dtrailer, rand(1:P.n))
        push!(ddriver, rand(1:P.k))
    elseif length(trailer) > 1
        pop!(dtrailer)
        pop!(ddriver)
    else
        return trailer, driver
    end
    s = gets(trailer, dtrailer, driver, ddriver, P)
    prob = kimura(s, P.N)
    if prob > rand()
        return dtrailer, ddriver

    else
        return trailer, driver

    end
end


"""
    gets(trailer1::Array{Int,1}, trailer2::Array{Int,1},
            driver1::Array{Int,1}, driver2::Array{Int,1}, P::parameters)

Compute selection coefficients after length mutation.
"""
function gets(trailer1::Array{Int,1}, trailer2::Array{Int,1},
            driver1::Array{Int,1}, driver2::Array{Int,1}, P::parameters)

    OldE = geten(trailer1, driver1)
    NewE = geten(trailer2, driver2)
    OldFit = P.f(OldE, length(trailer1), P)
    NewFit = P.f(NewE, length(trailer2), P)
    s = NewFit - OldFit
    return s
end


"""
    Step(trailer::Array{Int,1}, driver::Array{Int,1}, P::parameters)

Perform a step in the Gillespie Algorithm.
"""
function Step(trailer::Array{Int,1}, driver::Array{Int,1}, P::parameters)
    l = length(driver)
    r1 = P.N * (P.mu + P.nu) + P.rho
    r2 = rand()
    if r2 <= P.N*P.mu/r1
        trailer = tmut(trailer, driver, P)
    elseif r2 <= (P.N*P.mu+P.N*P.nu)/r1
        trailer, driver = lmut(trailer, driver, P)
    else
        driver = dmut(driver, P)
    end
    return trailer, driver
end


"""
    sim_steps(runtime::Integer, P::parameters)

Perform n runs of `runtime` steps.
"""
function sim_steps(runtime::Integer, n, P::parameters)
    E = zeros(n)
    L = zeros(n)
    for i in 1:n
        T = 0.
        t = 0.
        P.l = 10

        trailer, driver = createPop(P)
        while T <= runtime
            trailer, driver = Step(trailer, driver, P)
            T += 1
        end
        E[i], L[i] = geten(trailer,driver), length(trailer)
    end
    return E, L
end
