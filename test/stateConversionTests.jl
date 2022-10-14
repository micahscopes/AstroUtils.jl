
using AstroUtils, StaticArrays

# Test 1: Eliptical non-equatorial

# Position and velocity
rv = SVector(6524.834, 6862.875,  6448.296, # [km] 
             4.901327, 5.533756, -1.976341) # [km/s]

# Gravitational parameter
μ = 398600.4418 # [km³/s²]

# Convert
(kep,flag) = convertState(rv, Cartesian, Keplerian, μ)

# Tests 
@test abs(kep[1] - 36127.3376) < 0.01
@test abs(kep[2] - 0.8328534)  < 0.0001
@test abs(kep[3] - 1.5336056)  < 0.001
@test abs(kep[4] - 3.9775750)  < 0.001
@test abs(kep[5] - 0.9317428)  < 0.0001
@test abs(kep[6] - 1.6115525)  < 0.001
@test flag == 0
