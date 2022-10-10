
include("CubicSpline.jl")
include("Ephemeride.jl")

struct Ephemerides 
    # Initial and final epochs of ephemerides (in TDB)
    t0::Float64
    tf::Float64

    # Target IDs
    targIDs::Vector{Int}

    # Vector of Ephemeris data
    ephems::Vector{Ephemeride}
end

function Ephemerides(tspan, nPoints, targIDs, obsID, ref; abcorr = "NONE")
    # Create vector of Ephemerides
    ephems = [Ephemeride(tspan, nPoints, targIDs[i], obsID, ref; abcorr = abcorr) for i in eachindex(targIDs)]

    return Ephemerides(tspan[1], tspan[2], targIDs, ephems)
end

Base.eachindex(itr::Ephemerides) = eachindex(itr.ephems)
Base.getindex(r::Ephemerides, i) = r.ephems[i]

getTargetIDs(ephems::Ephemerides) = ephems.targIDs

function getGM(ephems::Ephemerides, targID)
    if !(targID in ephems.targIDs)
        throw(ArgumentError("Ephemerides does not contain ephemeris for the requested target ID."))
    end

    for i in eachindex(ephems)
        if getTargetID(ephems.ephems[i]) == targID
            return getTargetGM(ephems[i])
        end
    end
end

function getState(ephems::Ephemerides, targID, t)
    if !(targID in ephems.targIDs)
        throw(ArgumentError("Ephemerides does not contain ephemeris for the requested target ID."))
    end

    for i in eachindex(ephems)
        if getTargetID(ephems.ephems[i]) == targID
            return getState(ephems[i], t)
        end
    end
end

function getPosition(ephems::Ephemerides, targID, t)
    if !(targID in ephems.targIDs)
        throw(ArgumentError("Ephemerides does not contain ephemeris for the requested target ID."))
    end

    for i in eachindex(ephems)
        if getTargetID(ephems.ephems[i]) == targID
            return getPosition(ephems[i], t)
        end
    end
end

export Ephemeride, Ephemerides
export getTargetIDs, getGM, getState, getPosition