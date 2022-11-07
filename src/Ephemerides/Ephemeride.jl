
struct Ephemeride 
    # Initial and final epochs of ephemeride (in TDB)
    t0::Float64
    tf::Float64

    # Target NAIF ID
    targ::Int

    # Target gravitational parameter
    mu::Float64

    # Observer NAIF ID
    obs::Int

    # Refernece frame
    frame::String

    # Aberration corrections
    abcorr::String

    # Cubic spline interpolant
    spline::CubicSpline
end

function Ephemeride(tspan, nPoints, targ, obs, ref; abcorr = "NONE")
    # Compute epochs
    ts  = LinRange(0.0, 1.0, nPoints)

    # Get body gravitational parameter
    mu  = bodvcd(targ, "GM", 1)[1]

    # Get body's state at each epoch
    states = zeros(nPoints, 6)
    for i in 1:nPoints 
        state, ld   = spkez(targ, tspan[1] + (tspan[2] - tspan[1])*ts[i], ref, abcorr, obs)
        states[i,:] .= state
    end

    # Compute cubic spline interpolant
    spline = CubicSpline(ts, states)

    return Ephemeride(tspan[1], tspan[2], targ, mu, obs, ref, abcorr, spline)
end

# Simple getters for settings used to generate ephemeris
getTargetGM(ephem::Ephemeride)              = ephem.mu
getTargetID(ephem::Ephemeride)              = ephem.targ
getObserverID(ephem::Ephemeride)            = ephem.obs
getReferenceFrame(ephem::Ephemeride)        = ephem.frame
getAberrationCorrections(ephem::Ephemeride) = ephem.abcorr

function getState(ephem::Ephemeride, t)
    # Check that t is in [t0, tf]
    if t < ephem.t0 || t > ephem.tf
        throw(ArgumentError("Time passed to getState() which is outside of the ephemeride's coverage."))
    end

    # Unscale time
    ts = (t - ephem.t0)/(ephem.tf - ephem.t0)

    # Get state from spline and return
    return getState(ephem.spline, ts) 
end

function getPosition(ephem::Ephemeride, t)
    # Check that t is in [t0, tf]
    if t < ephem.t0 || t > ephem.tf
        throw(ArgumentError("Time passed to getState() which is outside of the ephemeride's coverage."))
    end

    # Unscale time
    ts = (t - ephem.t0)/(ephem.tf - ephem.t0)

    # Get state from spline and return
    return getPosition(ephem.spline, ts) 
end

function getPositionPartial(ephem::Ephemeride, t)
    # Check that t is in [t0, tf]
    if t < ephem.t0 || t > ephem.tf
        throw(ArgumentError("Time passed to getPositionPartial() which is outside of the ephemeride's coverage."))
    end

    # Unscale time
    ts = (t - ephem.t0)/(ephem.tf - ephem.t0)

    # Get state from spline and return
    dsdts = getPositionPartial(ephem.spline, ts)
    dtsdt = 1.0 / (ephem.tf - ephem.t0)
    return SVector(dsdts[1]*dtsdt, dsdts[2]*dtsdt, dsdts[3]*dtsdt)
end