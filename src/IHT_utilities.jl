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
    gk   :: Vector{Float64}     # Numerator of step size μ.   gk = df[idx] is a temporary array of length `k` that arises as part of the gradient calculations. I avoid doing full gradient calculations since most of `b` is zero. 
    xgk  :: Vector{Float64}     # Demonimator of step size μ. x * gk also part of the gradient calculation 
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
    b    = zeros(Float64, p)
    df   = zeros(Float64, p)
    xb   = zeros(Float64, n)
    r    = zeros(Float64, n)
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

"""
this function updates the BitArray indices for b. 
"""
function _iht_indices(
    v :: IHTVariable,
    k :: Int
)
    # set v.idx[i] = 1 if v.b[i] != 0 (i.e. find components of beta that are non-zero)
    v.idx .= v.b .!= 0

    # if idx is the 0 vector, v.idx[i] = 1 if i is one of the k largest components
    # of the gradient (stored in v.df), and set other components of idx to 0. 
    if sum(v.idx) == 0
        a = select(v.df, k, by=abs, rev=true) 
        v.idx[abs.(v.df) .>= abs(a)-2*eps()] .= true
        v.gk .= zeros(sum(v.idx))
    end

    return nothing
end

# """
#     fill_perm!(x, y, idx) 

# This subroutine fills a `k`-vector `x` from a `p`-vector `y` via an index vector `idx`.
# This variant admits BitArray index vectors.

# Arguments:

# - `x` is the `k`-vector to fill.
# - `y` is the `p`-vector to use in filling `x`.
# - `idx` is either a `BitArray` or `Int` vector` that indexes the components of `y` to put into `x`. If `idx` contains `Int`s, then only the first `k` indices are used. Otherwise, `fill_perm!()` traverses `idx` until it encounters `k` `true`s.
# """
# function fill_perm!(
#     x   :: Vector{Float64},
#     y   :: Vector{Float64},
#     idx :: BitArray{1}
# )
#     # x should have one element per "true" in idx
#     k = length(x)
#     #@assert k == sum(idx)
    
#     # counter j is used to track the number of trues in idx
#     j = 0

#     # loop over entire vector idx
#     @inbounds for i in eachindex(idx) 

#         # if current component of idx is a true, then increment j and fill x from y
#         if idx[i]
#             j += 1
#             x[j] = y[i]
#         end

#         # once x has k components, then it is completely filled and we return it
#         j == k && return nothing
#     end

#     return nothing
# end

"""
    project_k!(x, k)

This function projects a vector `x` onto the set S_k = { y in R^p : || y ||_0 <= k }.
It does so by first finding the pivot `a` of the `k` largest components of `x` in magnitude.
`project_k!` then thresholds `x` by `abs(a)`, sending small components to 0. 

Arguments:

- `b` is the vector to project.
- `k` is the number of components of `b` to preserve.
"""
function project_k!(
    x :: Vector{Float64},
    k :: Int;
)
    a = select(x, k, by = abs, rev = true)
    for i in eachindex(x) 
        if abs(x[i]) < abs(a) 
            x[i] = 0.0
        end
    end
    return nothing
end


function compute_ω(
    v         :: IHTVariable,
    snpmatrix :: Matrix{Float64}, 
    μ         :: Float64,
    k         :: Int
)
    # In order to compute ω, need β^{m+1} and xβ^{m+1}. Following eq.5,
    BLAS.axpy!(μ, v.df, v.b) # take the gradient step: v.b = β - μ∇f(β)
    project_k!(v.b, k)       # P_k( β - μ∇f(β) ): preserve top k components of b
    _iht_indices(v, k)       # Update idx. (find indices of new beta that are nonzero)
    sum(v.idx) <= k || warn("More than k components of b is non-zero! Need: VERY DANGEROUS DARK SIDE HACK!")

    #compute xβ^{m+1} based on β^{m+1} just calculated 
    A_mul_B!(v.xb, snpmatrix, v.b)

    #calculate ω efficiently (old b0 and xb0 have been copied before calling iht!)
    return sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)
end