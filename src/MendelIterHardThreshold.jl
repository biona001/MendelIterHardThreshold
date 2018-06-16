module MendelIterHardThreshold

using MendelBase
using SnpArrays
using IHT
using DataFrames

export IterHardThreshold

include("IHT_data_structure.jl")

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

    #
    #convert bitarrays to Float64 genotype matrix, and add a column of ones for intercept
    #
    snpmatrix = convert(Array{Float64,2}, x.snpmatrix)
    snpmatrix = [snpmatrix ones(size(snpmatrix, 1))]

    #
    # Begin IHT calculations
    #
    fill!(v.xb, 0.0) #initialize β = 0 vector, so Xβ = 0
    copy!(v.r, y)    #redisual = y-Xβ = y

    #compute the gradient once, and update v.xb to get the loop started
    BLAS.gemv!('T', -1.0, snpmatrix, v.r, 1.0, v.df) # save -X'(y - Xβ) to v.df

    for mm_iter = 1:max_iter
        
        # calculate the gradient ∇f(β)
        BLAS.gemv!('T', -1.0, snpmatrix, v.r, 1.0, v.df) # save -X'(y - Xβ) to v.df

        #calculate the step size μ 
        (mu, mu_step) = iht!(v, snpmatrix, y, k, nstep=max_step, iter=mm_iter)

        # save values from previous iterate
        copy!(v.b0, v.b)   # b0 = b
        copy!(v.xb0, v.xb) # Xb0 = Xb
        loss = next_loss

        return mu_step
    end


    return v.df
end #function L0_reg

"""
Calculates the IHT step size, and number of times line search was done. 
"""
function iht!(
    v     :: IHTVariable,
    x     :: Matrix{Float64},
    y     :: Vector{Float64},
    k     :: Int;
    iter  :: Int = 1,
    nstep :: Int = 50,
)
    μ = norm(v.b, 2)^2 / norm(v.df, 2)^2 
    ω = norm(v.b - v.b0, 2)^2 / norm(v.xb - v.xb0, 2)^2

    println(v.b)
    println(v.b0)

    println(v.xb)
    println(v.xb0)

    println(μ)
    println(ω)

    return (μ, ω)

    μ_step = 0
    for i = 1:nstep
        if μ < ω; break; end
        μ /= 2 #step halving

        # warn if mu falls below machine epsilon
        μ <= eps(typeof(μ)) && warn("Step size equals zero, algorithm may not converge correctly")

        μ_step += 1
    end

    return (μ, μ_step)
end


end # module

