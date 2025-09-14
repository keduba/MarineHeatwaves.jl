

"""
    a,b,r2,sigma_a,sigma_b,sigma_e = linreg(x,y)

Linear regression

y ≈ a + b *  x

"""
function linreg(x::AbstractVector, y::AbstractVector)

    # remove NaNs
    ind = .!isnan.(x) .& .!isnan.(y)
    x = x[ind]
    y = y[ind]

    xm = mean(x)
    xa = x .- xm
    ym = mean(y)
    ya = y .- ym

    ssxy = xa' * ya
    ssxx = xa' * xa
    ssyy = ya' * ya

    b = ssxy / ssxx
    a = ym - b * xm
    r2 = ssxy^2 / (ssxx * ssyy)

    n = length(x)

    sigma_e = sqrt(max((ssyy - b^2 * ssxx) / (n - 2), 0))
    sigma_a = sigma_e * sqrt(sum(x .^ 2) / (n * ssxx))
    sigma_b = sigma_e / sqrt(ssxx)


    # Calculate the F-statistic
    f_stat = (ssxy^2 / ssxx) / (ssyy - ssxy^2 / ssxx) * (n - 2)

    # Calculate the p-value for the F-statistic
    p_value = 1 - cdf(FDist(1, n - 2), f_stat)

    # Calculate the t-statistic for the slope (b)
    t_stat_b = b / sigma_b

    # Calculate the p-value for the t-statistic
    p_value_b = 2 * (1 - cdf(TDist(n - 2), abs(t_stat_b)))

    # Calculate T-Distribution p-value
    t_stat = b / sigma_b
    p_value_t = 2 * (1 - cdf(TDist(n - 2), abs(t_stat)))

    # Calculate F-Distribution p-value
    f_stat = r2 / (1 - r2) * (n - 2) / 1
    p_value_f = 1 - cdf(FDist(1, n - 2), f_stat)
    return a, b, r2, sigma_a, sigma_b, sigma_e
end

"""
Multi-dimensional linear regression
y ≈ a + X b

a is a scalar and b a vector

y = a + X b + ϵ

Σ ϵ² = (y- (a + X b))^2


y = X beta + ϵ

beta = inv(X*X') * X' * y
"""
function linreg(X::AbstractMatrix, y::AbstractVector)
    N = size(X, 1)
    n = size(X, 2)
    sX = sum(X, dims=1)
    sy = sum(y)
    # check for dimension match in X and y
    meanY = sy / N

    M = zeros(eltype(X), n + 1, n + 1)
    M[1, 1] = N
    M[1, 2:n+1] = sX
    M[2:n+1, 1] = sX
    M[2:n+1, 2:n+1] = X' * X

    beta = M \ vcat(sy, X' * y)

    a = beta[1]
    b = beta[2:end]

    # https://en.wikipedia.org/w/index.php?title=Fraction_of_variance_unexplained&oldid=913207135
    SSerr = zero(eltype(X))
    SStot = zero(eltype(X))
    for i = 1:N
        # prediction
        yp = a + (@view X[i, :])' * b

        SStot += (y[i] - meanY)^2
        SSerr += (y[i] - yp)^2
    end

    #@show SStot, SSerr
    # fraction of variance unexplained
    FVU = SSerr / SStot
    R2 = 1 - FVU
    return a, b, R2
end

