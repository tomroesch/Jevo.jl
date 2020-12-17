using Jevo, Test, BenchmarkTools, LinearAlgebra

@testset "Initiating " begin
    @testset "Simple binding sites" begin

        pop = Jevo.binding_sites(N=4, l=10, n=4)
        Jevo.initiate!(pop, 3)
        # Test population size
        @test sum(pop.freqs) == pop.N
        # Test number of species
        @test length(pop.seqs) == 3
        # Test sequence length of species
        @test length.(pop.seqs) == (ones(Int64, 3) .* pop.l)

        emat = ones(4,4)
        pop = Jevo.binding_sites(N=4, l=10, n=4)
        @test_throws DimensionMismatch Jevo.initiate!(pop, emat)

        emat = ones(4, 10)
        emat[1, :] .-= 1
        pop = Jevo.binding_sites(N=4, l=10, n=4)
        Jevo.initiate!(pop, emat)
        @test pop.seqs[1] == ones(Int, 10)



    end

    @testset "Driver Trailer sites" begin
        @testset "Random sequences" begin

            pop1 = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
            Jevo.initiate!(pop1, 3)
            # Test population size
            @test sum(pop1.freqs) == pop1.N
            # Test number of species
            @test length(pop1.seqs) == 3
            # Test sequence length of species
            @test length.(pop1.seqs) == (ones(Int64, 3) .* pop1.L)
            # Test length of driver
            @test length(pop1.driver) == pop1.L
            #@test_nowarn Jevo.initiate!(pop1, overwrite=true)

        end

        @testset "Optimal Sequence" begin
            pop = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
            emat = ones(4, 3)
            @test_throws DimensionMismatch Jevo.initiate!(pop, emat)


            pop = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
            emat = ones(4, 4)
            emat[1, :] .-= 1
            Jevo.initiate!(pop, emat)
            @test pop.seqs[1] == ones(Int, 10)
            @test length(pop.seqs) == 1

            emat = ones(4, 4) - Matrix{Float64}(I, 4, 4)
            pop = Jevo.driver_trailer_l(N=4, l_0=10, n=4, m=4, L=20)
            Jevo.initiate!(pop, emat)
            @test pop.seqs[1] == pop.driver
            @test length(pop.freqs) == 1
            @test pop.freqs[1] == 4

        end

        @testset "Given Driver" begin

            pop2 = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
            Jevo.initiate!(pop2, 3, driver=[1, 2, 3, 4])
            # Test that partial driver sequence is incorporated
            @test pop2.driver[1:4] == [1, 2, 3, 4]
            pop3 = Jevo.driver_trailer(N=4, L=10, n=4, m=12)
            Jevo.initiate!(pop3, 3, driver=collect(1:12))
            # Test that driver too long is cut short
            @test pop3.driver == collect(1:10)
            # Test warning message
            @test_logs (:warn, "l>L. Choosing L=l.") Jevo.driver_trailer(N=4, l=10, L=2, n=4)

            pop4 = Jevo.driver_trailer(N=4, L=10, n=4, m=4)
            # Test error message for wrong driver sequence
            @test_throws ArgumentError Jevo.initiate!(pop4, 3, driver=collect(1:12))

        end
    end

    @testset "Driver Trailer sites variable length" begin
        @testset "Random Sequences" begin

            pop1 = Jevo.driver_trailer_l(N=4, l_0=10, n=4, m=4, L=20)
            Jevo.initiate!(pop1, 3)
            # Test population size
            @test sum(pop1.freqs) == pop1.N
            # Test number of species
            @test length(pop1.seqs) == 3
            # Test sequence length of species
            @test length.(pop1.seqs) == (ones(Int64, 3) .* pop1.L)
            # Test length of driver
            @test length(pop1.driver) == pop1.L
            # test binding site lengths
            @test length(pop1.l) == 3
            @test pop1.l == ones(Int64, 3) * pop1.l_0

        end

        @testset "Optimal Sequence" begin
            pop = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
            emat = ones(4, 3)
            @test_throws DimensionMismatch Jevo.initiate!(pop, emat)

            emat = ones(4, 4)
            emat[1, :] .-= 1
            pop = Jevo.driver_trailer_l(N=4, l_0=10, n=4, m=4, L=20)
            Jevo.initiate!(pop, emat)
            @test pop.seqs[1] == ones(Int, 20)
            @test length(pop.seqs) == 1
            
            emat = ones(4, 4) - Matrix{Float64}(I, 4, 4)
            pop = Jevo.driver_trailer_l(N=4, l_0=10, n=4, m=4, L=20)
            Jevo.initiate!(pop, emat)
            @test pop.seqs[1] == pop.driver
            @test length(pop.freqs) == 1
            @test pop.freqs[1] == 4
            @test length(pop.l) == 1
            @test pop.l[1] == 10
        end

        @testset "Given Trailer" begin
            pop2 = Jevo.driver_trailer_l(N=4, l_0=10, n=4, m=4)
            Jevo.initiate!(pop2, 3, driver=[1, 2, 3, 4])
            # Test that partial driver sequence is incorporated
            @test pop2.driver[1:4] == [1, 2, 3, 4]

            pop3 = Jevo.driver_trailer_l(N=4, L=10, n=4, m=12)
            Jevo.initiate!(pop3, 3, driver=collect(1:12))
            # Test that driver too long is cut short
            @test pop3.driver == collect(1:10)
            # Test warning message
            @test_logs (:warn, "l_0>L. Choosing L=l_0.") Jevo.driver_trailer_l(N=4, l_0=10, L=2, n=4)

            pop4 = Jevo.driver_trailer_l(N=4, L=10, n=4, m=4)
            # Test error message for wrong driver sequence
            @test_throws ArgumentError Jevo.initiate!(pop4, 3, driver=collect(1:12))
        end
    end

    @testset "Monomorphic Population" begin

        @testset "Initiate Random" begin
            pop = Jevo.mono_pop(N=10, l=10, n=4, m=4)
            Jevo.initiate!(pop)
            @test length(pop.seqs) == 10
            @test length(pop.driver) == 10
        end

        @testset "Initiate Optimal" begin
            pop = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
            emat = ones(4, 3)
            @test_throws DimensionMismatch Jevo.initiate!(pop, emat)

            emat = ones(4, 4)
            emat[1, :] .-= 1
            pop = Jevo.mono_pop(N=10, l=10, n=4, m=4)
            Jevo.initiate!(pop, emat)
            @test length(pop.seqs) == 10
            @test length(pop.driver) == 10
            @test pop.seqs == ones(Int64, 10)

            
            emat = ones(4, 4) - Matrix{Float64}(I, 4, 4)
            pop = Jevo.mono_pop(N=10, l=10, n=4, m=4)
            Jevo.initiate!(pop, emat)
            @test length(pop.seqs) == 10
            @test length(pop.driver) == 10
            @test pop.seqs == pop.driver
        end

    end
end


@testset "Mutations" begin
    @testset "Find existing sequences" begin
        pop = Jevo.binding_sites(N=8, l=1, n=4, seqs=[[1], [2], [3], [4]], freqs=[2, 2, 2, 2])
        Jevo.mutation!(pop)
        # Test that mutation to existing species does not add new species
        @test length(pop.seqs) == 4
        @test length(pop.freqs) == 4
    end
    @testset "Add new species" begin
        pop = Jevo.binding_sites(N=2, l=1, n=2)
        Jevo.initiate!(pop, 1)
        Jevo.mutation!(pop)
        # Test that new species is created
        @test length(pop.seqs) == 2
        @test length(pop.freqs) == 2
    end
    @testset "Driver mutation" begin
        pop = Jevo.driver_trailer(N=4, L=10, n=4, m=12)
        Jevo.initiate!(pop, 1, driver=collect(1:12))
        Jevo.mutation!(pop)
        # Test that new species is created
        @test length(pop.seqs) == 2
        @test length(pop.freqs) == 2

        pop2 = Jevo.driver_trailer(N=4, L=10, n=4, m=4)
        Jevo.initiate!(pop2, 3)
        tmp_driver = copy(pop2.driver)
        Jevo.driver_mutation!(pop2)
        # Test that driver sequence is mutated
        @test tmp_driver != pop2.driver
        @test sum(tmp_driver .!= pop2.driver) == 1
    end
    @testset "Length mutation" begin
        pop = Jevo.driver_trailer_l(N=4, L=10, n=4, m=12)
        Jevo.initiate!(pop, 1, driver=collect(1:12))
        Jevo.mutation!(pop)
        # Test that new species is created
        @test length(pop.seqs) == 2
        @test length(pop.freqs) == 2
        @test length(pop.l) == 2

        pop2 = Jevo.driver_trailer_l(N=4, L=10, n=4, m=4)
        Jevo.initiate!(pop2, 3)
        tmp_driver = copy(pop2.driver)
        Jevo.driver_mutation!(pop2)
        # Test that driver sequence is mutated
        @test tmp_driver != pop2.driver
        @test sum(tmp_driver .!= pop2.driver) == 1

        pop3 = Jevo.driver_trailer_l(N=1, L=10, n=4, m=4)
        Jevo.initiate!(pop3, 1)
        Jevo.mutation!(pop3)
        # Test that new species is created
        @test length(pop3.seqs) == 1
        @test length(pop3.freqs) == 1
        @test length(pop3.l) == 1


        pop4 = Jevo.driver_trailer_l(N=1, L=10, n=4, m=4)
        Jevo.initiate!(pop4, 1)
        Jevo.length_mutation!(pop4)
        # Test that new species is created
        @test length(pop4.seqs) == 1
        @test length(pop4.freqs) == 1
        @test length(pop4.l) == 1
    end
end

@testset "Substitutions" begin
    pop = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
    Jevo.initiate!(pop, 2)
    emat = ones(4, 4) .- Matrix{Float64}(I, 4, 4)
    f = Jevo.fermi_fitness(l=10, beta=1, f0=1, fl=0)
    # Test error for too many species
    @test_throws ArgumentError Jevo.bp_substitution!(pop, emat, f)

    pop = Jevo.driver_trailer(N=4, l=10, n=4, m=4)
    Jevo.initiate!(pop, 1)
    temp_seqs = deepcopy(pop.seqs)

    Jevo.bp_substitution!(pop, emat, f)

    @test length(pop.seqs) == 1
    @test length(pop.seqs[1]) == 10

    @test sum(temp_seqs[1] .!= pop.seqs[1]) <= 1

    @testset "Monomorphic Population" begin

    end
end
