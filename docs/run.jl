using SnpArrays
x = SnpArray("test")
snpmatrix = convert(Array{Float64,2}, x)

for i in 1:size(snpmatrix, 2)
	snpmatrix[:, i] = (snpmatrix[:, i] .- mean(snpmatrix[:, i])) / std(snpmatrix[:, i])
end


using SnpArrays
srand(2018)
x = SnpArray(rand(0:2, 5, 3))
view(x, :, [true, false, true])
#how to do the following?
view(x, :, [true, false, true]) * rand(2)




using IHT
# x = MendelIHT("test_control.txt")
x = MendelIHT("gwas 1 Control.txt")

using MendelIterHardThreshold
# x = IterHardThreshold("test_control.txt")
x = IterHardThreshold("gwas 1 Control.txt")




function plus()
	s = 0
	@simd for i = 1:10^6
		@inbounds s += i
	end
	return s
end

function plus_new()
	s = 0
	for i = 1:10^6
		@inbounds s += i
	end
	return s
end

#benchmark std computation
using SnpArrays, BenchmarkTools
x = SnpData("gwas 1 data")
std_vec = zeros(x.snps)
function hi()
    for i in 1:x.snps
        storage = view(x.snpmatrix, :, i)
        std_vec[i] .= 1.0 / std( convert(Vector{Float64}, storage) )
    end
	return std_vec
end 

storage = zeros(x.people)
function hii()
    for i in 1:x.snps
    	copy!(storage, x.snpmatrix[:, i]) #in place version of convert
        std_vec[i] .= 1.0 ./ std(storage)
    end
	return std_vec
end 

@benchmark hi()
@benchmark hii()


function hiii()
    mean_vec = 2.0x.maf #multiply by 2 because mean of each snp = 2.0*maf
    std_vec = 1.0 ./ sqrt.(2.0 .* x.maf .* (1.0 .- x.maf))
    return mean_vec, std_vec
end



using SnpArrays, BenchmarkTools
x = SnpArray("gwas 1 data")
z = bitrand(10000)
y = SnpArray(2200, sum(z))
function hii()
	x = SnpArray("gwas 1 data")
	z = bitrand(10000)
	y = SnpArray(2200, sum(z))
	y .= view(x, :, z)
end


using SnpArrays, PLINK
x = SnpData("test")
y = BEDFile("test.bed")
x.snpmatrix
y[:, 1:2] #obmit 3rd column, which is the intercept
convert(Matrix{Float64}, x.snpmatrix)
