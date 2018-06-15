module MendelIterHardThreshold

using MendelBase
using SnpArrays
using IHT
using DataFrames

export IterHardThreshold

"""
This is the wrapper function for the Iterative Hard Thresholding analysis option. 
"""
function IterHardThreshold(control_file = ""; args...)
    const IHT_VERSION :: VersionNumber = v"0.1.0"
    #
    # Print the logo. Store the initial directory.
    #
    print(" \n \n")
    println("     Welcome to OpenMendel's")
    println("      IHT analysis option")
    println("        version ", IHT_VERSION)
    print(" \n \n")
    println("Reading the data.\n")
    initial_directory = pwd()
    #
    # The user specifies the analysis to perform via a set of keywords.
    # Start the keywords at their default values.
    #
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    #
    # Define some keywords unique to this analysis option. 
    #
    keyword["data_type"] = ""
    keyword["predictors"] = ""
    #
    # Process the run-time user-specified keywords that will control the analysis.
    # This will also initialize the random number generator.
    #
    process_keywords!(keyword, control_file, args)
    #
    # Check that the correct analysis option was specified.
    #
    lc_analysis_option = lowercase(keyword["analysis_option"])
    if (lc_analysis_option != "" &&
      lc_analysis_option != "iht")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
    end
    keyword["analysis_option"] = "Iterative Hard Thresholding"
    #
    # Read the genetic data from the external files named in the keywords.
    #
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
    #
    # Execute the specified analysis.
    #
    println(" \nAnalyzing the data.\n")
##
    result = iht_gwas(person, snpdata, pedigree_frame, keyword)
    return result
##
    # execution_error = iht_gwas(person, snpdata, pedigree_frame, keyword)
    # if execution_error
    #   println(" \n \nERROR: Mendel terminated prematurely!\n")
    # else
    #   println(" \n \nMendel's analysis is finished.\n")
    # end
    # #
    # # Finish up by closing, and thus flushing, any output files.
    # # Return to the initial directory.
    # #
    # close(keyword["output_unit"])
    # cd(initial_directory)
    # return nothing
end #function IterHardThreshold

"""
This function performs IHT on GWAS data. The overall strategy: Make IHTVariables, use that to 
rewrite L0_reg, and modify iht!. Need to rewrite 3 functions in 3 weeks. 
"""
function iht_gwas(person::Person, snpdata::SnpData,
  pedigree_frame::DataFrame, keyword::Dict{AbstractString, Any})
    phenotype = convert(Array{Float64,1}, pedigree_frame[:Trait])
    
    x = L0_reg(snpdata, phenotype, 10)
    return x
end

"""
Object to contain intermediate variables and temporary arrays. Used for cleaner code in L0_reg
"""
type IHTVariable
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

"""
Regression step
"""
function L0_reg(
    x        :: SnpData,
    y        :: Vector{Float64}, 
    k        :: Int;
    v        :: IHTVariable = IHTVariables(x, y, k),
    mask_n   :: Vector{Int} = ones(Int, size(y)),
    tol      :: Float64 = 1e-4,
    max_iter :: Int     = 100,
    max_step :: Int     = 50,
)
    # first handle errors
    k        >= 0            || throw(ArgumentError("Value of k must be nonnegative!\n"))
    max_iter >= 0            || throw(ArgumentError("Value of max_iter must be nonnegative!\n"))
    max_step >= 0            || throw(ArgumentError("Value of max_step must be nonnegative!\n"))
    tol      >  eps(Float64) || throw(ArgumentError("Value of global tol must exceed machine precision!\n"))
    
    # initialize return values
    mm_iter   = 0                 # number of iterations of L0_reg
    mm_time   = 0.0               # compute time *within* L0_reg
    next_loss = oftype(tol,Inf)   # loss function value

    # initialize floats
    current_obj = oftype(tol,Inf) # tracks previous objective function value
    the_norm    = 0.0             # norm(b - b0)
    scaled_norm = 0.0             # the_norm / (norm(b0) + 1)
    mu          = 0.0             # Landweber step size, 0 < tau < 2/rho_max^2

    # initialize integers
    i       = 0                   # used for iterations in loops
    mu_step = 0                   # counts number of backtracking steps for mu

    # initialize booleans
    converged = false             # scaled_norm < tol?

    #convert bitarrays to Float64 genotype matrix, and add a column of ones for intercept
    snpmatrix = convert(Array{Float64,2}, x.snpmatrix)
    snpmatrix = [snpmatrix ones(size(snpmatrix, 1))]

    # update xb, r, and gradient
    # initialize_xb_r_grad!(v, x, y, k, pids=pids)
    if sum(v.idx) == 0
        fill!(v.xb, 0.0)
        copy!(v.r, y)
        v.r[mask_n .== 0] = 0.0
    else
        #A_mul_B!(v.xb, x, v.b, v.idx, k, mask_n, pids=pids)
        A_mul_B!(v.xb, x, v.b, v.idx, k, mask_n)
        difference!(v.r, y, v.xb)
        v.r[mask_n .== 0] = 0.0
    end

    for mm_iter = 1:max_iter
        
        # calculate the gradient
        BLAS.gemv!('T', -1.0, snpmatrix, v.r, 1.0, v.df) # -X'(y - Xβ)

        (mu, mu_step_halving) = iht!(v, snpmatrix, y, k, nstep=max_step, iter=mm_iter)

        println(mu)
        return mu
    end


    return v.df
end #function L0_reg

"""
The IHT step
"""
function iht!(
    v     :: IHTVariable,
    x     :: Matrix{Float64},
    y     :: Vector{Float64},
    k     :: Int;
    iter  :: Int = 1,
    nstep :: Int = 50,
)
  return (12, 13)
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


end # module

