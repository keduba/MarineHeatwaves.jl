using Random

@testset "det mylabel" begin
    vec_data = collect(Float32, repeat(1.:5., 2))
    vec_exceeds = MarineHeatwaves._excess(vec_data, 3)
    function create_test(data::VecOrMat, exceeds::BitArray)
        MHW(
            data,  # temp
            zeros(Float32, size(data)),  # clim (placeholder)
            zeros(Float32, size(data)),  # thresh (placeholder)
            exceeds,
            EventHW
        )
    end
    ms = create_test(vec_data, vec_exceeds)

    @test_throws ErrorException MarineHeatwaves.mylabel(ms, Int16(5), Int16(2))
    @test_throws ErrorException MarineHeatwaves.mylabel(ms, Int16(-1), Int16(2))
end

@testset "det anomsa" begin
    # Some data
    ms = rand(Xoshiro(1), Float32, 20, 5)
    mc = rand(Xoshiro(2), Float32, 20, 5)
    mh = rand(Xoshiro(3), Float32, 20, 5)
    me = MHW(ms, mc, mh)
    starts = Vector{Int16}[
        Int16[1, 4, 10, 15],
        Int16[2, 5, 8, 15]
    ]
    ends = Vector{Int16}[
        Int16[3, 5, 13, 20],
        Int16[4, 8, 12, 20]
    ]
    cols = [1, 2]
    enst = (starts, ends, cols)
    ahoy = MarineHeatwaves.anomsam(me, enst)
    typeof(ahoy) == MarineHeatwaves.EventsFull
    length(ahoy.means) == length(cols)
    size(ahoy.tpanom) == (size(ms, 1), length(cols))
end
