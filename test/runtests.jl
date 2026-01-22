using MarineHeatwaves
using Test
using Dates

@testset "MarineHeatwaves.jl" begin
    # Write your tests here.
    println("Testing ...")
    @test leapyearday(Date(2013, 3,31)) == 91
end

ad = Date(2012, 11, 3): Date(2014, 8, 26)
cr = Date(2013, 9, 3): Date(2014, 1, 20)
jk = MarineHeatwaves.timeindices(ad, cr)

@testset "timefunctions.jl" begin
    # Write your tests here.
    println("Testing time functions...")
    ad = Date(2012, 11, 3): Date(2014, 8, 26)
    cr = Date(2013, 9, 3): Date(2014, 1, 20)
    jk = MarineHeatwaves.timeindices(ad, cr)
    
    @test leapyearday(Date(2013, 3,31)) == 91
    @test jk[1] == 305:444
    @test length(jk[1]) == length(jk[2]) == 140
    @test length(MarineHeatwaves.daterange(jk[2], 3)) == 366
end


@testset "mhwarrays.jl" begin
    println("Testing the array functions... ")
	tmp = rand(50,50,50)
	@test MarineHeatwaves.seamask(rand(50), 1) == BitArray(1)
	@test ndims(MarineHeatwaves._subtemp(tmp, 1:30)[1])  == 2

	
	# drg = MarineHeatwaves.daterange(jk[2], 3)
    # @test MarineHeatwaves.tresh(indt, drg, 0.5)
	#@test size(MarineHeatwaves.clim(rand(5000), drg), 1) == 366
    # @test MarineHeatwaves.exceed()
end

@testset "Full Tests" begin
    @testset "" begin
    	
    end
end
