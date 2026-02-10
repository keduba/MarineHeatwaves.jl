using Test
using Dates

# Test structure (simple)
#  * correct working
#  * error-handling
# test leap years
# test date/day ranges

@testset "deterministic - time functions" begin
    # Proper functioning
    println("Testing time functions...")
    main_start = Date(2012, 6, 3)
    main_end = Date(2012, 11, 26)
    sub_start = Date(2012, 9, 3)
    sub_end = Date(2012, 9, 20)
    ab = main_start:main_end
    cr = sub_start:sub_end
    cd = main_start:sub_start
    de = sub_end:main_end
    jk = MarineHeatwaves.timeindices(ab, cr)
    @test leapyearday(Date(2013, 3, 31)) == 91
    @test jk[1] == 93
    @test length(jk) == length(cr) == 18
    @test jk isa UnitRange

    # Error handling
    @test_throws MethodError leapyearday("stringy")
    # the main event should be longer
    @test_throws ErrorException MarineHeatwaves.timeindices(cr, ab)
    @test_throws ErrorException MarineHeatwaves.timeindices(cd, de)
end

@testset "property-based -- time functions" begin
end
