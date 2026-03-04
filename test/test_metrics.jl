
function create_test_events(n_events::Int, dtype::Type)
    # Helper function to create a test EventsFull struct
    T = Float32
    Events(
        [T.(rand(1:100, n_events)) for _ in 1:1],   # means
        [T.(rand(1:100, n_events)) for _ in 1:1],   # minimaxes
        [T.(rand(1:100, n_events)) for _ in 1:1],   # onset
        [T.(rand(1:100, n_events)) for _ in 1:1],   # decline
        [T.(rand(1:10, n_events)) for _ in 1:1],    # duration
        [T.(rand(1:1000, n_events)) for _ in 1:1],  # sums
        [T.(rand(1:5, n_events)) for _ in 1:1],   # categorys
        rand(T, n_events, 1),             # tpanom
        rand(T, n_events, 1),             # thanom
        dtype
    )
end

function create_test_events(n_pixels::Int, n_events_per_pixel::Int, dtype::Type)
    # Helper function to create a test EventsFull struct
    T = Float32
    Events(
        [T.(rand(1:100, n_events_per_pixel)) for _ in 1:n_pixels],   # means
        [T.(rand(1:100, n_events_per_pixel)) for _ in 1:n_pixels],   # minimaxes
        [T.(rand(1:100, n_events_per_pixel)) for _ in 1:n_pixels],   # onset
        [T.(rand(1:100, n_events_per_pixel)) for _ in 1:n_pixels],   # decline
        [T.(rand(1:10, n_events_per_pixel)) for _ in 1:n_pixels],    # duration
        [T.(rand(1:1000, n_events_per_pixel)) for _ in 1:n_pixels],  # sums
        [T.(rand(1:5, n_events_per_pixel)) for _ in 1:n_pixels],   # categorys
        rand(T, n_events_per_pixel, n_pixels),             # tpanom
        rand(T, n_events_per_pixel, n_pixels),             # thanom
        dtype
    )
end

metrics = MarineHeatwaves.metrics

@testset "det meanmetrics" begin
    dtp = EventHW
    emv = create_test_events(5, dtp)
    mdate = [Date(2020, 1, 1), Date(2021, 1, 1), Date(2022, 1, 1), Date(2023, 1, 1), Date(2024, 1, 1)]
    result = MarineHeatwaves.meanmetricsm(emv, mdate)
    @test size(result, 1) == length(getfield(emv, :means))
    @test size(result, 2) == length(metrics)
    # Verify frequency calculation
    @test result[metrics[:frequency]] ≈ (length.(emv.duration) ./ 5)[]  atol=1e-6
    # Verify days calculation
    @test result[metrics[:days]] ≈ (sum.(emv.duration) / 5)[] atol=1e-6
end

@testset "det annualmetrics" begin
    dtp = EventHW
    n_px = 3
    n_events_per_px = 5
    emv = create_test_events(n_px, n_events_per_px, dtp)
    # Create a step range of dates
    mdate = Date(2020,1,1):Month(1):Date(2021,12,31)

    evst = (cst=[1:2], cse=[3:4], cols=[1,2])

    result = MarineHeatwaves.annualmetricsm(emv, mdate, evst)

    @test any(!isnan, result)
    @test any(result[:,:,metrics[:frequency]] .> 0)
    @test any(result[:,:,metrics[:days]] .> 0)
end

@testset "det annualmetrics single year" begin
        # Create test data for a single year
        dtype = EventCS
        n_pixels = 2
        n_events_per_pixel = 4
        emv = create_test_events(n_pixels, n_events_per_pixel, dtype)
        
        # Create a step range for a single year
        mdate = Date(2020,1,1):Month(1):Date(2020,12,31)
        
        # Event start type
        evst = (cst=[1:2], cse=[3:4], cols=[1,2])
        
        # Call the function
        result = MarineHeatwaves.annualmetricsm(emv, mdate, evst)
        
        # Verify result dimensions for single year
        @test size(result, 1) == 1  # single year
    @test size(result, 2) == length(evst[3])  # pixels/columns
    @test size(result, 3) == length(metrics)  # metrics
end

@testset "trend functions" begin
    n_years = 5; n_pixels = 3; n_metrics = length(metrics)
    outannual = zeros(Float32, n_years, n_pixels, n_metrics)
    for m in 1:n_metrics
	for p in 1:n_pixels
	    base = rand() * 10
	    slope = rand() * 0.5 - 0.25
	    for y in 1:n_years
	        outannual[y, p, m] = base + slope * y + randn() * 0.1
	    end
	end
    end
    outannual[:, 1, metrics[:means]] = [1.0, 2.0, 3.0, 4.0, 5.0]
    outannual[:, 2, metrics[:sums]] = [5.0, 4.0, 3.0, 2.0, 1.0]
    
    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue = trendm(outannual)
    @test outcoeff[1, metrics[:sums]] > 0
    @test all(0 .<= outrsqd .<= 1)
    @test all(0 .<= outpvalue .<= 1)
end

@testset "special trend" begin
    n_years = 5; n_pixels = 3; n_metrics = length(metrics)
    outannual = zeros(Float32, n_years, n_pixels, n_metrics)

    for p in 1:n_pixels
	outannual[:, p, metrics[:means]] = [1.0, 2.0, 3.0, 4.0, 5.0] .+ p
    end

    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue = trendm(outannual)
    @test all(outcoeff[:, metrics[:means]] .> 0)
    @test all(outrsqd[:, metrics[:means]] .> 0.9)
end
