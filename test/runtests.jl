using AstroUtils
using Test, SafeTestsets

@time begin
@time @safetestset "State conversion tests..." begin include("stateConversionTests.jl") end
@time @safetestset "State conversion partials tests..." begin include("stateConversionPartialsTests.jl") end
end