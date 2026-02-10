# using Test

@testset "deterministic -- range" begin
    # start code here
    @test begin
        TI = Int16
        clyd = collect(TI, 1:10)
        win::TI = 3
        ranges = MarineHeatwaves.urange(clyd, win)
        ranges[1] == [1:4]
        length(ranges) == 366
    end
end

@testset "property-based -- range" begin
end

@testset "deterministic -- seamask" begin
    @test begin
        # "1-D Array"
        test_1d = [1., 2., NaN, 4.]
        seamask(test_1d, 1)[1] == 1
    end
    
    @test begin
        # "3-D Array"
        test_3d = rand(7, 7, 7)
        test_3d[3:5, 3:5, 3:5] .= NaN
        size(seamask(test_3d, 3)) == (7, 7)
    end
end

@testset "property-based -- seamask" begin
end

@testset "deterministic -- subtemp" begin
    # Testing 1-D array
    sst_1d = collect(Float32, 1:10)
    mhix = 2:6
    clix = 1:7

    res = _subtemp(sst_1d, mhix, clix)
    @test length(res) == 3
    @test res[1] == sst_1d[mhix]
    @test res[2] == sst_1d[clix]

    # Testing 3-D Array
    sst_3d = rand(5, 5, 6)
    mhix = 1:3
    clix = 2:5
    ms, cs, cix = _subtemp(sst_3d, mhix, clix)
    @test length(cix) == 4
    @test size(ms) == (length(mhix), length(cix[1]))
    @test size(cs) == (length(clix), length(cix[1]))
end

@testset "deterministic -- moving_mean" begin
    vec = collect(1.:8.)

    TI = Int16
    # 3-point moving mean, no wrap: wrap = false
    mm = MarineHeatwaves.moving_mean(vec, TI(3), false)
    @test mm ≈ [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.5]
    # 3-point moving mean, wrap = true
    mm = MarineHeatwaves.moving_mean(vec, TI(3), true)
    @test mm[2:7] ≈ [2.0, 3.0, 4.0, 5., 6.0, 7.0]
end

@testset "deterministic -- moving_means" begin
    vec = collect(1.:8.)
    TI = Int16
    ly = 1:8
    # 3-point moving mean, wrap = true
    mm = MarineHeatwaves.moving_mean(vec, TI(3), true)
    @test mm[2:7] ≈ [2.0, 3.0, 4.0, 5., 6.0, 7.0]
    mms = MarineHeatwaves.moving_means(vec, TI(3), ly, wrap=true)
    @test mm == mms
    @test length(mms) == length(ly)
end

@testset "det -- climthresh" begin
    println("Testing climthresh")
    ngrid = 3
    cmst = Matrix{Float32}(reshape(1.0:366*ngrid, 366, ngrid))
    mly = Int16[60, 150, 250] #1:200
    cly = collect(Int16, 1:366)
    win = 5
    swin = 31
    thv = 0.9f0

    cr, tr = MarineHeatwaves.climthresh(cmst, cly, mly, win, swin, thv, wrap=true)
    @test size(cr) == (length(mly), ngrid)
    @test size(tr) == (length(mly), ngrid)
end

