using MendelIterHardThreshold
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
include("MendelIterHardThreshold_test.jl")
include("IHT_utilities_test.jl")

# julia -e 'Pkg.test("MendelIterHardThreshold",coverage=true)'
# @show get_summary(process_file("src/MendelIterHardThreshold.jl"))
# @show get_summary(process_file("src/IHT_utilities.jl"))