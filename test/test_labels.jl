
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
    # data first
    # md = rand(8, 8)
    # ahoy = anomsa(m, enst, ixs)
    # typeof(ahoy) == EventFull
end
