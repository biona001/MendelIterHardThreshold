module MendelIterHardThreshold

using MendelBase
using SnpArrays
using IHT
using DataFrames
using Distances
using StatsBase

export IterHardThreshold

include("IHT_utilities.jl")

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
    
    k = keyword["predictors"]

    x = L0_reg(snpdata, phenotype, k)
    return x
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

    # start timer
    tic()

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

    #convert bitarrays to Float64 genotype matrix, standardize each SNP, and add intercept
    snpmatrix = convert(Array{Float64,2}, x.snpmatrix)
    # snpmatrix = use_A2_as_minor_allele(x.snpmatrix) #to compare results with using PLINK
    snpmatrix = StatsBase.zscore(snpmatrix, 1) 
    snpmatrix = [ones(size(snpmatrix, 1)) snpmatrix] 

    #
    # Begin IHT calculations
    #
    fill!(v.xb, 0.0) #initialize β = 0 vector, so Xβ = 0
    copy!(v.r, y)    #redisual = y-Xβ = y
    v.r[mask_n .== 0] .= 0 #bit masking? idk why we need this yet

    # calculate the gradient v.df = X'(y - Xβ) one time. Future gradient calculations are done in iht!
    # Can we use v.xk instead of snpmatrix?
    # Why aren't we doing v.df = -X'(y - Xβ)?
    At_mul_B!(v.df, snpmatrix, v.r) 

    for mm_iter = 1:max_iter
        # save values from previous iterate
        copy!(v.b0, v.b)   # b0 = b
        copy!(v.xb0, v.xb) # Xb0 = Xb
        loss = next_loss
        
        #calculate the step size μ. Can we use v.xk instead of snpmatrix?
        (μ, μ_step) = iht!(v, snpmatrix, y, k, nstep=max_step, iter=mm_iter)

        # iht! gives us an updated x*b. Use it to recompute residuals and gradient
        v.r .= y .- v.xb
        v.r[mask_n .== 0] .= 0 #bit masking, idk why we need this yet 

        At_mul_B!(v.df, snpmatrix, v.r) # v.df = X'(y - Xβ) Can we use v.xk instead of snpmatrix?

        # update loss, objective, gradient, and check objective is not NaN or Inf
        next_loss = sum(abs2, v.r) / 2
        !isnan(next_loss) || throw(error("Objective function is NaN, aborting..."))
        !isinf(next_loss) || throw(error("Objective function is Inf, aborting..."))

        # track convergence
        the_norm    = chebyshev(v.b, v.b0) #max(abs(x - y))
        scaled_norm = the_norm / (norm(v.b0, Inf) + 1)
        converged   = scaled_norm < tol

        if converged
            # stop time
            mm_time = toq()

            println("IHT converged in " * string(mm_iter) * " iterations")
            println("It took " * string(mm_time) * " seconds to converge")
            println("The estimated model is: " * string(v.b))
            println("There are " * string(countnz(v.b)) * " non-zero entries of β")
            return v.b
        end

        if mm_iter == max_iter
            # stop time
            mm_time = toq()

            println("Did not converge!!!!! omg!!!!!")
            println("It took " * string(mm_time) * " seconds to run IHT (in vain)")

            return nothing
        end
    end
end #function L0_reg

"""
Calculates the IHT step β+ = P_k(β - μ ∇f(β)). 
Returns step size (μ), and number of times line search was done (μ_step). 

This function updates: b, xb, xk, gk, xgk, idx
"""
function iht!(
    v         :: IHTVariable,
    snpmatrix :: Matrix{Float64},
    y         :: Vector{Float64},
    k         :: Int;
    iter      :: Int = 1,
    nstep     :: Int = 50,
)
    # compute indices of nonzeroes in beta and store them in v.idx (also sets size of v.gk)
    _iht_indices(v, k)

    # fill v.xk, which stores columns of snpmatrix corresponding to non-0's of b
    v.xk[:, :] .= snpmatrix[:, v.idx]

    # fill v.gk, which store only k largest components of gradient (v.df)
    # fill_perm!(v.gk, v.df, v.idx)  # gk = g[v.idx]
    v.gk .= v.df[v.idx]

    # now compute X_k β_k and store result in v.xgk
    A_mul_B!(v.xgk, v.xk, v.gk)

    # warn if xgk only contains zeros
    all(v.xgk .≈ 0.0) && warn("Entire active set has values equal to 0")

    #compute step size and notify if step size too small
    μ = norm(v.gk, 2)^2 / norm(v.xgk, 2)^2 
    isfinite(μ) || throw(error("Step size is not finite, is active set all zero?"))
    μ <= eps(typeof(μ)) && warn("Step size $(μ) is below machine precision, algorithm may not converge correctly")

    #Take the gradient step and compute ω. Note in order to compute ω, need β^{m+1} and xβ^{m+1} (eq5)
    _iht_gradstep(v, μ, k)
    ω = compute_ω!(v, snpmatrix) #is snpmatrix required? Or can I just use v.x

    #compute ω and check if μ < ω. If not, do line search by halving μ and checking again.
    μ_step = 0
    for i = 1:nstep
        #exit loop if μ < ω where c = 0.01 for now
        if _iht_backtrack(v, ω, μ); break; end 

        #if μ >= ω, step half and warn if μ falls below machine epsilon
        μ /= 2 
        μ <= eps(typeof(μ)) && warn("Step size equals zero, algorithm may not converge correctly")

        # recompute gradient step
        copy!(v.b, v.b0)
        _iht_gradstep(v, μ, k)

        # re-compute ω based on xβ^{m+1}
        A_mul_B!(v.xb, snpmatrix, v.b)
        ω = sqeuclidean(v.b, v.b0) / sqeuclidean(v.xb, v.xb0)

        μ_step += 1
    end

    return (μ, μ_step)
end


end # module

