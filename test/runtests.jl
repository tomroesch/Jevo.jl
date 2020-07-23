using Jedi, Test, BenchmarkTools

@testset "Initiating " begin
    @testset "Simple binding sites" begin
        pop = Jedi.binding_sites(N=4, l=10, n=4)
        Jedi.initiate_rand!(pop, 3)
        # Test population size
        @test sum(pop.freqs) == pop.N
        # Test number of species
        @test length(pop.seqs) == 3
        # Test sequence length of species
        @test length.(pop.seqs) == (ones(Int64, 3) .* pop.L)
        # Test warning message
        @test_logs (:warn, "l>L. Choosing L=l.") Jedi.binding_sites(N=4, l=10, L=2, n=4)
    end

    @testset "Driver Trailer sites" begin
        pop1 = Jedi.driver_trailer(N=4, l=10, n=4, m=4)
        Jedi.initiate_rand!(pop1, 3)
        # Test population size
        @test sum(pop1.freqs) == pop1.N
        # Test number of species
        @test length(pop1.seqs) == 3
        # Test sequence length of species
        @test length.(pop1.seqs) == (ones(Int64, 3) .* pop1.L)
        # Test length of driver
        @test length(pop1.driver) == pop1.L

        pop2 = Jedi.driver_trailer(N=4, l=10, n=4, m=4)
        Jedi.initiate_rand!(pop2, 3, driver=[1, 2, 3, 4])
        # Test that partial driver sequence is incorporated
        @test pop2.driver[1:4] == [1, 2, 3, 4]

        pop3 = Jedi.driver_trailer(N=4, L=10, n=4, m=12)
        Jedi.initiate_rand!(pop3, 3, driver=collect(1:12))
        # Test that driver too long is cut short
        @test pop3.driver == collect(1:10)
        # Test warning message
        @test_logs (:warn, "l>L. Choosing L=l.") Jedi.driver_trailer(N=4, l=10, L=2, n=4)

        pop4 = Jedi.driver_trailer(N=4, L=10, n=4, m=4)
        # Test error message for wrong driver sequence
        @test_throws ArgumentError Jedi.initiate_rand!(pop4, 3, driver=collect(1:12))
    end
end


@testset "Mutations" begin
    @testset "Find existing sequences" begin
        pop = Jedi.binding_sites(N=8, L=1, n=4, seqs=[[1], [2], [3], [4]], freqs=[2, 2, 2, 2])
        Jedi.mutation!(pop)
        # Test that mutation to existing species does not add new species
        @test length(pop.seqs) == 4
        @test length(pop.freqs) == 4
    end
    @testset "Add new species" begin
        pop = Jedi.binding_sites(N=2, L=1, n=2)
        Jedi.initiate_rand!(pop, 1)
        Jedi.mutation!(pop)
        # Test that new species is created
        @test length(pop.seqs) == 2
        @test length(pop.freqs) == 2
    end
    @testset "Driver mutation" begin
        pop = Jedi.driver_trailer(N=4, L=10, n=4, m=12)
        Jedi.initiate_rand!(pop, 1, driver=collect(1:12))
        Jedi.mutation!(pop)
        # Test that new species is created
        @test length(pop.seqs) == 2
        @test length(pop.freqs) == 2

        pop2 = Jedi.driver_trailer(N=4, L=10, n=4, m=4)
        Jedi.initiate_rand!(pop2, 3)
        tmp_driver = copy(pop2.driver)
        Jedi.driver_mutation!(pop2)
        # Test that driver sequence is mutated
        @test tmp_driver != pop2.driver
        @test sum(tmp_driver .!= pop2.driver) == 1
    end
end
