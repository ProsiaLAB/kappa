module Opacities
include("Mie.jl")
include("Components.jl")
include("Fractal.jl")
include("Geofractal.jl")

include("Utils.jl")

using .Components
using .Fractal
using .Geofractal
using .Mie
using .Utils

g = Mie.add_two(3)
println(g)
printstyled(g, color=:red)

export Components, Fractal, Geofractal, Mie, Utils

end # module Opacities
