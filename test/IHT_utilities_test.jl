using MendelIterHardThreshold, SnpArrays, MendelBase, CSV, DataFrames

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

function test_data()
	#dataset with 2 SNP and 6 people. The SNP matrix is 6x3 (with column of intercept)
	x = SnpData("test") 
	y = CSV.read("test.fam", delim = ' ', header = false)
	y = convert(Array{Float64,1}, y[:, 6])
	k = 2
	v = IHTVariables(x, y, 2)
	return (x, y, k, v)
end

function gwas1_data()
	#dataset with 10000 SNP and 2200 people. The SNP matrix is 2200x10001 (with column of intercept)
	x = SnpData("test") 
	y = CSV.read("test.fam", delim = ' ', header = false)
	y = convert(Array{Float64,1}, y[:, 6])
	k = 10
	v = IHTVariables(x, y, 2)
	return (x, y, k, v)
end


@testset "initilize IHTVariables" begin
	(x, y, k, v) = test_data()

	#k must be an integer, but could be larger than 
	@test_throws(ArgumentError, IHTVariables(x, y, 0))
	@test_throws(ArgumentError, IHTVariables(x, y, -1))
	@test_throws(MethodError, IHTVariables(x, y, 1.1))
	@test_throws(MethodError, IHTVariables(x, y, NaN))
	@test_throws(MethodError, IHTVariables(x, y, missing))
	@test_throws(MethodError, IHTVariables(x, y, Inf))

	#Different types of inputs for IHTVariables(x, y, k) is 
	@test typeof(v) == IHTVariable
	@test typeof(x) == SnpData || typeof(x) <: SnpArray

	@test size(v.b)    == (3,) 
	@test size(v.b0)   == (3,)
	@test size(v.xb)   == (6,)
	@test size(v.xb0)  == (6,)
	@test size(v.xk)   == (6, 2)
	@test size(v.gk)   == (2,)
	@test size(v.xgk)  == (6,)
	@test size(v.idx)  == (3,)
	@test size(v.idx0) == (3,)
	@test size(v.r)	   == (6,)
	@test size(v.df)   == (3,)

	@test typeof(v.b)    == Array{Float64, 1}
	@test typeof(v.b0)   == Array{Float64, 1}
	@test typeof(v.xb)   == Array{Float64, 1}
	@test typeof(v.xb0)  == Array{Float64, 1}
	@test typeof(v.xk)   == Array{Float64, 2}
	@test typeof(v.gk)   == Array{Float64, 1}
	@test typeof(v.xgk)  == Array{Float64, 1}
	@test typeof(v.idx)  == BitArray{1}
	@test typeof(v.idx0) == BitArray{1}
	@test typeof(v.r)	 == Array{Float64, 1}
	@test typeof(v.df)   == Array{Float64, 1}
end

@testset "_iht_indices" begin
	(x, y, k, v) = gwas1_data()
end

@testset "project_k!" begin
    
end

@testset "compute_ω!" begin
    
end

@testset "use_A2_as_minor_allele" begin
    
end

@testset "_iht_backtrack" begin
    
end

@testset "_iht_gradstep" begin
    
end