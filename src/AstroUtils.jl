module AstroUtils

using LinearAlgebra
using StaticArrays
using SparseArrays
using SPICE
using IfElse
using Downloads: download

const webLSK    = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const webSPK    = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const webPCK    = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/Gravity.tpc"
const pkgSrcDir     = @__DIR__
const defaultLSK    = pkgSrcDir * "/../data/kernels/lsk/naif0012.tls"
const defaultSPK    = pkgSrcDir * "/../data/kernels/spk/de440.bsp"
const defaultPCK    = pkgSrcDir * "/../data/kernels/pck/gravity.tpc"

# Download default spice kernels if necessary
!isfile(defaultLSK) ? download(webLSK, defaultLSK) : ()
!isfile(defaultSPK) ? download(webSPK, defaultSPK) : ()
!isfile(defaultPCK) ? download(webPCK, defaultPCK) : ()

# Includes
include("spice.jl")
include("StateRepresentations/stateRepresentations.jl")
include("Ephemerides/Ephemerides.jl")

export furnshDefaults

end
