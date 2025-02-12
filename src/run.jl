using ArgParse
#using Profile
#using PProf

#--------------------------------------------------------
# The problem name is a command line argument:
#
# 1. Launch Julia:
# >> julia --project=.
#
# 2. Push equations name to ARGS
#    You need this only when you run a new equations
#
#    julia > push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);
#    julia > include(./src/Jexpresso.jl)
#
#    EQUATIONS is the name of your equations directory as $JEXPRESSO/src/equations/EQUATIONS
#    EQUATIONS_CASE_NAME is the name of the subdirectory $JEXPRESSO/src/equations/EQUATIONS_CASE_NAME
#
# Ex. To run the Compressible Euler equations in $JEXPRESSO/src/equations/CompEuler/theta
# 
#  julia > push!(empty!(ARGS), "CompEuler", "theta");
#  julia > include(./src/Jexpresso.jl)
#
#--------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "eqs"
        help = "equations"
        default = "CompEuler"
        required = false
        
        "eqs_case"
        help = "case name in equations directory"
        default = "theta"
        required = false
        
    end

    return parse_args(s)
end

#--------------------------------------------------------
#Parse command line args:
#--------------------------------------------------------
parsed_args                = parse_commandline()
parsed_equations           = string(parsed_args["eqs"])
parsed_equations_case_name = string(parsed_args["eqs_case"])

driver_file          = string(dirname(@__DIR__()), "/problems/equations/drivers.jl")
case_name_dir        = string(dirname(@__DIR__()), "/problems/equations", "/", parsed_equations, "/", parsed_equations_case_name)
user_input_file      = string(case_name_dir, "/user_inputs.jl")
user_flux_file       = string(case_name_dir, "/user_flux.jl")
user_source_file     = string(case_name_dir, "/user_source.jl")
user_bc_file         = string(case_name_dir, "/user_bc.jl")
user_initialize_file = string(case_name_dir, "/initialize.jl")

include(driver_file)

include(user_input_file)
include(user_flux_file)
include(user_source_file)
include(user_bc_file)
include(user_initialize_file)

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
mod_inputs_print_welcome()
inputs = Dict{}()

inputs = mod_inputs_user_inputs!(user_input_file)

#--------------------------------------------------------
#Create output directory if it doesn't exist:
#--------------------------------------------------------
user_defined_output_dir = inputs[:output_dir]
if user_defined_output_dir == "none"
    OUTPUT_DIR = string(case_name_dir, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
    inputs[:output_dir] = OUTPUT_DIR
else
    OUTPUT_DIR = string(user_defined_output_dir, "/", parsed_equations, "/", parsed_equations_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
end
if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

#--------------------------------------------------------
#Save a copy of user_inputs.jl for the case being run 
#--------------------------------------------------------
run(`$cp $user_input_file $OUTPUT_DIR`) 


driver(ContGal(),   # Space discretization type    
       inputs, # input parameters from src/user_input.jl
       OUTPUT_DIR,
       TFloat)

# Export pprof profile and open interactive profiling web interface.
#pprof()
