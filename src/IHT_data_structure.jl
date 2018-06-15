export IHTVariable, IHTVariables

"""
Object to contain intermediate variables and temporary arrays. Used for cleaner code in L0_reg
"""
struct IHTVariable
    b    :: Vector{Float64}     # the statistical model, most will be 0
    b0   :: Vector{Float64}     # previous estimated model in the mm step
    xb   :: Vector{Float64}     # vector that holds x*b 
    xb0  :: Vector{Float64}     # previous xb in the mm step
    xk   :: Matrix{Float64}     # the n by k subset of the design matrix x corresponding to non-0 elements of b
    gk   :: Vector{Float64}     # gk = df[idx] is a temporary array of length `k` that arises as part of the gradient calculations. I avoid doing full gradient calculations since most of `b` is zero.
    xgk  :: Vector{Float64}     # x * gk also part of the gradient calculation
    idx  :: BitArray{1}         # BitArray indices of nonzeroes in b for A_mul_B
    idx0 :: BitArray{1}         # previous iterate of idx
    r    :: Vector{Float64}     # n-vector of residuals
    df   :: Vector{Float64}     # the gradient: df = -x' * (y - xb)

    # IHTVariable(b::DenseVector, b0::Vector{Float64}, xb::DenseVector, 
    #   xb0::Vector{Float64}, xk::Matrix{Float64}, gk::Vector{Float64}, 
    #   xgk::Vector{Float64}, idx::BitArray{1}, idx0::BitArray{1}, r::DenseVector, 
    #   df::DenseVector) = new{Float64,DenseVector}(b, b0, xb, xb0, xk, gk, xgk, idx, idx0, r, df)
end

function IHTVariables(
    x :: SnpData,
    y :: Vector{Float64},
    k :: Int64
)
    n    = x.people
    p    = x.snps + 1 #adding 1 because we need an intercept
    # pids = procs(x)
    # V    = typeof(y)
    b    = Vector{Float64}(p,) 
    df   = Vector{Float64}(p,) 
    xb   = Vector{Float64}(n,) 
    r    = Vector{Float64}(n,) 
    b0   = zeros(Float64, p)
    xb0  = zeros(Float64, n)
    xk   = zeros(Float64, n, k)
    xgk  = zeros(Float64, n)
    gk   = zeros(Float64, k)
    idx  = falses(p) 
    idx0 = falses(p)
    return IHTVariable(b, b0, xb, xb0, xk, gk, xgk, idx, idx0, r, df)
end

#for testing purposes, keeping for future reference
#taken from read_plink_data in data.jl in PLINK
# function read_fam_file()
#   #
#   # read the FAM file
#     #
#     famfile = keyword["plink_input_basename"] * ".fam"
#     Y = readdlm(famfile, ' ', header=false)
#     #
#     # check that the FAM file has six column
#     #
#     p = size(Y,2)
#     p == 6 || throw(DimensionMismatch("FAM file does not have six columns, is it formatted correctly?"))
#     #
#     # in FAM file, the phenotype is the rightmost column 
#     # initialize a SharedVector and fill it with the phenotype.
#     # we cannot know for certain that Y is loaded as floating point
#     # to be safe, explicitly convert the phenotype column to the correct type
#     #
#     y = convert(Vector, Y[:,end])
#     return(y) 
# end