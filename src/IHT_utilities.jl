export IHTVariable, IHTVariables, use_A2_as_minor_allele, standardize_snpmatrix
export use_A2_as_minor_allele, project_k!, project_group_sparse!, _iht_gradstep, _init_iht_indices
export _iht_omega, _iht_backtrack, std_reciprocal, gIHTResults

"""
Object to contain intermediate variables and temporary arrays. Used for cleaner code in L0_reg
"""
mutable struct IHTVariable{T <: Float, V <: DenseVector}

    #TODO: Consider changing b and b0 to SparseVector
    #TODO: Add in itc
    # itc  :: T             # the intercept
    # itc0 :: T             # old intercept
    b     :: Vector{T}     # the statistical model, most will be 0
    b0    :: Vector{T}     # previous estimated model in the mm step
    xb    :: Vector{T}     # vector that holds x*b 
    xb0   :: Vector{T}     # previous xb in the mm step
    xk    :: SnpLike{2}    # the n by k subset of the design matrix x corresponding to non-0 elements of b
    gk    :: Vector{T}     # gk = df[idx] is a temporary array of length `k` that arises as part of the gradient calculations. I avoid doing full gradient calculations since most of `b` is zero. 
    xgk   :: Vector{T}     # x * gk also part of the gradient calculation 
    idx   :: BitVector     # BitArray indices of nonzeroes in b for A_mul_B
    idx0  :: BitVector     # previous iterate of idx
    r     :: V             # n-vector of residuals
    df    :: V             # the gradient: df = -x' * (y - xb - intercept)
    group :: Vector{Int64} # vector denoting group membership
end

function IHTVariables{T <: Float}(
    x :: SnpLike{2},
    y :: Vector{T},
    J :: Int64,
    k :: Int64
) 
    n, p = size(x) 
    p += 1 # add 1 for the intercept, need to change this to use itc later

    # itc  = zero(T)
    # itc0 = zero(T)
    b     = zeros(T, p)
    b0    = zeros(T, p)
    xb    = zeros(T, n)
    xb0   = zeros(T, n)
    xk    = SnpArray(n, J*k)
    gk    = zeros(T, J*k)
    xgk   = zeros(T, n)
    idx   = falses(p) 
    idx0  = falses(p)
    r     = zeros(T, n)
    df    = zeros(T, p)
    group = ones(Int64, p)

    return IHTVariable{T, typeof(y)}(b, b0, xb, xb0, xk, gk, xgk, idx, idx0, r, df, group)
    # return IHTVariable{T, typeof(y)}(itc, itc0, b, b0, xb, xb0, xk, gk, xgk, idx, idx0, r, df)
end

#"""
#Returns ω, a constant we need to bound the step size μ to guarantee convergence. 
#"""
#
#function compute_ω!(v::IHTVariable, snpmatrix::Matrix{Float64})
#    #update v.xb
#    A_mul_B!(v.xb, snpmatrix, v.b)
#
#    #calculate ω efficiently (old b0 and xb0 have been copied before calling iht!)
#    return sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)
#end

"""
This function is needed for testing purposes only. 
Converts a SnpArray to a matrix of float64 using A2 as the minor allele. We want this function 
because SnpArrays.jl uses the less frequent allele in each SNP as the minor allele, while PLINK.jl 
always uses A2 as the minor allele, and it's nice if we could cross-compare the results. 
"""
function use_A2_as_minor_allele(snpmatrix :: SnpArray)
    n, p = size(snpmatrix)
    matrix = zeros(n, p)
    for i in 1:p
        for j in 1:n
            if snpmatrix[j, i] == (0, 0); matrix[j, i] = 0.0; end
            if snpmatrix[j, i] == (0, 1); matrix[j, i] = 1.0; end
            if snpmatrix[j, i] == (1, 1); matrix[j, i] = 2.0; end
            if snpmatrix[j, i] == (1, 0); matrix[j, i] = missing; end
        end
    end
    return matrix
end

"""
This function computes the gradient step v.b = P_k(β + μ∇f(β)), and updates v.idx. 
Recall calling axpy! implies v.b = v.b + μ*v.df, but v.df stores an extra negative sign.
"""
function _iht_gradstep{T <: Float}(
    v  :: IHTVariable{T},
    μ  :: T,
    J  :: Int,
    k  :: Int;
)
   # n = length(v.xb)
   # v.itc = v.itc0 + μ*(n*μ - sum(v.r) + n*v.itc0)
   BLAS.axpy!(μ, v.df, v.b) # take the gradient step: v.b = b + μ∇f(b) (which is an addition since df stores X(-1*(Y-Xb)))
   project_group_sparse!(v.b, v.group, J, k) # project to doubly sparse vector 
   v.idx .= v.b .!= 0       # Update idx. (find indices of new beta that are nonzero)

   # If the k'th largest component is not unique, warn the user. 
   sum(v.idx) <= J*k || warn("More than J*k components of b is non-zero! Need: VERY DANGEROUS DARK SIDE HACK!")
end

"""
When initializing the IHT algorithm, take largest elements of each group of df as nonzero 
components of b. This function set v.idx = 1 for those indices. 
"""
function _init_iht_indices{T <: Float}(
    v     :: IHTVariable{T},
    J     :: Int,
    k     :: Int
)    
    perm = sortperm(v.df, by = abs, rev = true) 
    project_group_sparse!(v.df, v.group, J, k)
    v.idx[find(v.df)] = true
    v.gk = zeros(T, sum(v.idx))

    return nothing
end

"""
this function calculates the omega (here a / b) used for determining backtracking
"""
function _iht_omega{T <: Float}(
    v :: IHTVariable{T}
)
    a = sqeuclidean(v.b, v.b0::Vector{T}) :: T
    b = sqeuclidean(v.xb, v.xb0::Vector{T}) :: T
    return a, b
end

"""
this function for determining whether or not to backtrack. True = backtrack
"""
function _iht_backtrack{T <: Float}(
    v       :: IHTVariable{T},
    ot      :: T,
    ob      :: T,
    mu      :: T,
    mu_step :: Int,
    nstep   :: Int
)
    mu*ob > 0.99*ot              &&
    sum(v.idx) != 0              &&
    sum(xor.(v.idx,v.idx0)) != 0 &&
    mu_step < nstep
end

"""
Compute the standard deviation of a SnpArray in place
"""
function std_reciprocal{T <: Float}(A::SnpArray, mean_vec::Vector{T})
    m, n = size(A)
    @assert n == length(mean_vec) "number of columns of snpmatrix doesn't agree with length of mean vector"
    std_vector = zeros(T, n) 

    @inbounds for j in 1:n
        @simd for i in 1:m
            (a1, a2) = A[i, j]
            if !isnan(a1, a2) #only add if data not missing
                std_vector[j] += (convert(T, a1 + a2) - mean_vec[j])^2
            end 
        end
        std_vector[j] = 1.0 / sqrt(std_vector[j] / (m - 1))
    end
    return std_vector
end

""" Projects the point y onto the set with at most m active groups and at most 
n active predictors per group. The variable group encodes group membership. Currently
assumes there are no unknown group membership.
TODO: 1)make this function operate in-place 
      2)check if sortperm can be replaced by something that doesn't sort the whole array
"""
function project_group_sparse!{T <: Float}(
    y     :: Vector{T}, 
    group :: Vector{Int64}, 
    m     :: Int64, 
    n     :: Int64
)
    groups = maximum(group)
    group_count = zeros(Int, groups)         #counts number of predictors in each group
    z = zeros(groups)                        #l2 norm of each group
    perm = zeros(Int64, length(y))           #vector holding the permuation vector after sorting
    sortperm!(perm, y, by = abs, rev = true)

    #calculate the magnitude of each group, where only top predictors contribute
    for i in eachindex(y)
        j = perm[i] 
        k = group[j]
        if group_count[k] < n
            z[k] = z[k] + y[j]^2
            group_count[k] = group_count[k] + 1
        end
    end

    #go through the top predictors in order. Set predictor to 0 if criteria not met
    group_rank = zeros(Int64, length(z)) 
    sortperm!(group_rank, z, rev = true)
    group_rank = invperm(group_rank)
    fill!(group_count, 1) 
    for i in eachindex(y)
        j = perm[i]
        k = group[j]
        if (group_rank[k] > m) || (group_count[k] > n)
            y[j] = 0.0
        else
            group_count[k] = group_count[k] + 1
        end
    end
end

"""
an object that houses results returned from a group IHT run
"""
immutable gIHTResults{T <: Float, V <: DenseVector}
    time  :: T
    loss  :: T
    iter  :: Int
    beta  :: V
    J     :: Int64
    k     :: Int64
    group :: Vector{Int64}

    #gIHTResults{T,V}(time::T, loss::T, iter::Int, beta::V) where {T <: Float, V <: DenseVector{T}} = new{T,V}(time, loss, iter, beta)
    gIHTResults{T,V}(time, loss, iter, beta, J, k, group) where {T <: Float, V <: DenseVector{T}} = new{T,V}(time, loss, iter, beta, J, k, group)
end

# strongly typed external constructor for gIHTResults
gIHTResults(time::T, loss::T, iter::Int, beta::V, J::Int, k::Int, group::Vector{Int}) where {T <: Float, V <: DenseVector{T}} = gIHTResults{T, V}(time, loss, iter, beta, J, k, group)

"""
a function to display gIHTResults object
"""
function Base.show(io::IO, x::gIHTResults)
    println(io, "IHT results:")
    println(io, "\nCompute time (sec):     ", x.time)
    println(io, "Final loss:             ", x.loss)
    println(io, "Iterations:             ", x.iter)
    println(io, "Max number of groups:   ", x.J)
    println(io, "Max predictors/group:   ", x.k)
    println(io, "IHT estimated ", countnz(x.beta), " nonzero coefficients.")
    non_zero = find(x.beta)
    print(io, DataFrame(Group=x.group[non_zero], Predictor=non_zero, Estimated_β=x.beta[non_zero]))
    return nothing
end
