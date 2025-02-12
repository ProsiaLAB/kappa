module Opacities
include("Mie.jl")
include("Components.jl")
include("Fractal.jl")
include("Geofractal.jl")

include("Utils.jl")

using TOML

using .Components
using .Fractal
using .Geofractal
using .Mie
using .Utils

struct Params
    amin::Float64
    amax::Float64
    apow::Float64
    amean::Float64
    asig::Float64
    na::Int
    sdkind::String
    lmin::Float64
    lmax::Float64
    nlam::Int
    nang::Int
    chopangle::Float64
    pcore::Float64
    pmantle::Float64
    nm::Int
    nmant::Int
    method::String
    fmax::Float64
    write_fits::Bool
    write_scatter::Bool
    for_radmc::Bool
end



"""
    parse_config(file::String)

Loads and processes the given TOML configuration file.
"""
function parse_config(file::String)
    try
        config = TOML.parsefile(file)
        println("Successfully loaded configuration from $file")
        return config

        # Call processing functions here
        # process_config(config)

    catch e
        println("Error: Could not load file $file")
        println(e)
    end
end

"""
    calculate()

Automatically finds a TOML file in the current working directory and processes it.
"""
function calculate(file::String)

    if isnothing(file)
        println("No config file found in the current directory.")
        return
    end

    config = parse_config(file)

    # Main program
    mat_rho = zeros(Float64, 21)

    init_params(config)
end

"""
    init_params()

Initializes the parameters for the calculation and sets defaults.
"""
function init_params(config::Dict{String, Any})
    println(length(config["materials"]["entries"]))
end

end # module Opacities
