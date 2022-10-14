using AstroUtils, StaticArrays, ForwardDiff

# Position and velocity
rv = SVector(6524.834 + 100.0*randn(), 6862.875 + 100.0*randn(),  6448.296 + 100.0*randn(), # [km] 
             4.901327 + 1.0*randn(), 5.533756 + 1.0*randn(), -1.976341 + 1.0*randn()) # [km/s]
μ  = 3.986e5

# Convert to MEE 
mee = convertState(rv, Cartesian, MEE, μ)

# Compute partials with AstroUtils
cartWrtMee = convertStatePartials(mee, MEE, Cartesian, μ)

# Compute partials with ForwardDiff
cartWrtMeeFD = ForwardDiff.jacobian(x -> convertState(x, MEE, Cartesian, μ), mee)

# Compute diff
jacDiff = cartWrtMee .- cartWrtMeeFD
for i in eachindex(jacDiff)
    @test jacDiff[i] ./ eps(cartWrtMee[i]) < 10.0
end