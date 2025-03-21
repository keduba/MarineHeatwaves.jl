using Statistics
using Test

include("linreg.jl")


x = 2005:1/12:2022
y = 1 .+ 2*x


a,b,R2 = linreg(x,y)

@test a ≈ 1
@test b ≈ 2
@test R2 ≈ 1

y_pred = a .+ x*b

# use linreg with arrays (lon,lat,time)

a = randn(10,11)
b = randn(10,11)

Y = zeros(10,11,length(x))
for j = 1:size(Y,2)
    for i = 1:size(Y,1)
        Y[i,j,:] = a[i,j] .+ b[i,j]*x
    end
end


result = mapslices(slice -> linreg(x,slice),Y,dims=3)
result_a = getindex.(result,1)
result_b = getindex.(result,2)
result_R2 = getindex.(result,3)

@test result_a ≈ a
@test result_b ≈ b
@test all(result_R2 .≈ 1)

