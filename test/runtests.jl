using Jedi, Test, BenchmarkTools

@testset "Initiating " begin
    @testset "Simple binding sites" begin
        pop = Jedi.binding_sites(N=4, l=10, n=4)
        Jedi.initiate_rand!(pop, 3)
        @test sum(pop.freqs) == pop.N
        @test length(pop.seqs) == 3
        @test length.(pop.seqs) == (ones(Int64, 3) .* pop.L)
        @test_logs (:warn, "l>L. Choosing L=l.") Jedi.binding_sites(N=4, l=10, L=2, n=4)
    end

    @testset "Driver Trailer sites" begin
        pop1 = Jedi.driver_trailer(N=4, l=10, n=4, m=4)
        Jedi.initiate_rand!(pop1, 3)
        @test sum(pop1.freqs) == pop1.N
        @test length(pop1.seqs) == 3
        @test length.(pop1.seqs) == (ones(Int64, 3) .* pop1.L)
        @test length(pop1.driver) == pop1.L

        pop2 = Jedi.driver_trailer(N=4, l=10, n=4, m=4)
        Jedi.initiate_rand!(pop2, 3, driver=[1, 2, 3, 4])
        @test pop2.driver[1:4] == [1, 2, 3, 4]

        pop3 = Jedi.driver_trailer(N=4, L=10, n=4, m=4)
        Jedi.initiate_rand!(pop3, 3, driver=collect(1:12))
        @test pop3.driver == collect(1:10)

        @test_logs (:warn, "l>L. Choosing L=l.") Jedi.driver_trailer(N=4, l=10, L=2, n=4)
    end
end


@testset "Mutations" begin
    @testset "Find existing sequences" begin
        pop = Jedi.binding_sites(N=8, L=1, n=4, seqs=[[1], [2], [3], [4]], freqs=[2, 2, 2, 2])
        Jedi.mutation!(pop)
        @test length(pop.seqs) == 4
        @test length(pop.freqs) == 4
    end
    @testset "Add new species" begin
        pop = Jedi.binding_sites(N=2, L=1, n=2)
        Jedi.initiate_rand!(pop, 1)
        Jedi.mutation!(pop)
        @test length(pop.seqs) == 2
        @test length(pop.freqs) == 2
    end
end
