
convertStatePartials(state::AbstractArray, ::Type{MEE}, ::Type{Cartesian}, mu) = cartWrtMee(state, mu)
function cartWrtMee(mee::AbstractArray, mu)
    # Compute requirements
    ss      = 1.0 + mee[4]*mee[4] + mee[5]*mee[5]
    aa      = mee[4]*mee[4] - mee[5]*mee[5]
    w       = 1.0 + mee[2]*cos(mee[6]) + mee[3]*sin(mee[6])

    invP    = 1.0 / mee[1]
    invSs   = 1.0 / ss
    invW    = 1.0 / w
    onePaa  = 1.0 + aa
    oneMaa  = 1.0 - aa
    sqrmup  = sqrt(mu / mee[1])

    # ===== Position partials
    drxdp   = invSs*invW*(onePaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    drydp   = invSs*invW*(2.0*mee[4]*mee[5]*cos(mee[6]) + oneMaa*sin(mee[6]))
    drzdp   = 2.0*invSs*invW*(-mee[5]*cos(mee[6]) + mee[4]*sin(mee[6])) 

    drxdf   = -mee[1]*cos(mee[6])*invSs*invW*invW*(onePaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    drydf   = -mee[1]*cos(mee[6])*invSs*invW*invW*(2.0*mee[4]*mee[5]*cos(mee[6]) + oneMaa*sin(mee[6]))
    drzdf   = 2.0*mee[1]*cos(mee[6])*invSs*invW*invW*(mee[5]*cos(mee[6]) - mee[4]*sin(mee[6]))

    drxdg   = -mee[1]*sin(mee[6])*invSs*invW*invW*(onePaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    drydg   = mee[1]*sin(mee[6])*invSs*invW*invW*(-2.0*mee[4]*mee[5]*cos(mee[6]) - oneMaa*sin(mee[6]))
    drzdg   = 2.0*mee[1]*sin(mee[6])*invSs*invW*invW*(mee[5]*cos(mee[6]) - mee[4]*sin(mee[6]))

    drxdh   = 2.0*mee[5]*mee[1]*invSs*invSs*invW*(2.0*mee[4]*mee[5]*cos(mee[6]) + oneMaa*sin(mee[6]))
    drydh   = 2.0*mee[1]*invSs*invSs*invW*(mee[5]*oneMaa*cos(mee[6]) - 2.0*mee[4]*(1.0 + mee[5]*mee[5])*sin(mee[6]))
    drzdh   = 2.0*mee[1]*invSs*invSs*invW*(2.0*mee[4]*mee[5]*cos(mee[6]) + oneMaa*sin(mee[6]))

    drxdk   = 2.0*mee[1]*invSs*invSs*invW*(-2.0*(1.0 + mee[4]*mee[4])*mee[5]*cos(mee[6]) + mee[4]*onePaa*sin(mee[6]))
    drydk   = 2.0*mee[4]*mee[1]*invSs*invSs*invW*(onePaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    drzdk   = -2.0*mee[1]*invSs*invSs*invW*(onePaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))

    drxdL   = -mee[1]*invSs*invW*invW*(mee[3] + mee[3]*mee[4]*mee[4] - 2.0*mee[2]*mee[4]*mee[5] - mee[3]*mee[5]*mee[5] -
                2.0*mee[4]*mee[5]*cos(mee[6]) + onePaa*sin(mee[6]))
    drydL   = mee[1]*invSs*invW*invW*(mee[2] - mee[2]*mee[4]*mee[4] - 2.0*mee[3]*mee[4]*mee[5] + mee[2]*mee[5]*mee[5] + 
                oneMaa*cos(mee[6]) - 2.0*mee[4]*mee[5]*sin(mee[6]))
    drzdL   = 2.0*mee[1]*invSs*invW*invW*(mee[2]*mee[4] + mee[3]*mee[5] + mee[4]*cos(mee[6]) + mee[5]*sin(mee[6]))

    # ===== Velocity partials
    dvxdp   = 0.5*sqrmup*invSs*invP*(mee[3] + mee[3]*mee[4]*mee[4] - 2.0*mee[2]*mee[4]*mee[5] - mee[3]*mee[5]*mee[5] - 
                2.0*mee[4]*mee[5]*cos(mee[6]) + onePaa*sin(mee[6]))
    dvydp   = 0.5*sqrmup*invSs*invP*(-mee[2] + mee[2]*mee[4]*mee[4] + 2.0*mee[3]*mee[4]*mee[5] - mee[2]*mee[5]*mee[5] - 
                oneMaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    dvzdp   = -sqrmup*invSs*invP*(mee[2]*mee[4] + mee[3]*mee[5] + mee[4]*cos(mee[6]) + mee[5]*sin(mee[6]))

    dvxdf   = 2.0*mee[4]*mee[5]*sqrmup*invSs 
    dvydf   = oneMaa*sqrmup*invSs 
    dvzdf   = 2.0*mee[4]*sqrmup*invSs

    dvxdg   = -onePaa*sqrmup*invSs 
    dvydg   = -2.0*mee[4]*mee[5]*sqrmup*invSs 
    dvzdg   = 2.0*mee[5]*sqrmup*invSs

    dvxdh   = 2.0*mee[5]*sqrmup*invSs*invSs*(mee[2] - mee[2]*mee[4]*mee[4] - 2.0*mee[3]*mee[4]*mee[5] + mee[2]*mee[5]*mee[5] + 
                oneMaa*cos(mee[6]) - 2.0*mee[4]*mee[5]*sin(mee[6]))
    dvydh   = -2.0*sqrmup*invSs*invSs*(2.0*mee[2]*mee[4] + mee[3]*mee[5] - mee[3]*mee[4]*mee[4]*mee[5] + 2.0*mee[2]*mee[4]*mee[5]*mee[5] + 
                mee[3]*mee[5]*mee[5]*mee[5] + 2.0*mee[4]*(1.0 + mee[5]*mee[5])*cos(mee[6]) + mee[5]*oneMaa*sin(mee[6]))
    dvzdh   = -2.0*sqrmup*invSs*invSs*(-mee[2] + mee[2]*mee[4]*mee[4] + 2.0*mee[3]*mee[4]*mee[5] - mee[2]*mee[5]*mee[5] - 
                oneMaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    
    dvxdk   = 2.0*sqrmup*invSs*invSs*(mee[2]*mee[4] + mee[2]*mee[4]*mee[4]*mee[4] + 2.0*mee[3]*mee[5] + 2.0*mee[3]*mee[4]*mee[4]*mee[5] -
                mee[2]*mee[4]*mee[5]*mee[5] + mee[4]*onePaa*cos(mee[6]) + 2.0*(1.0 + mee[4]*mee[4])*mee[5]*sin(mee[6]))
    dvydk   = -2.0*mee[4]*sqrmup*invSs*invSs*(mee[3] + mee[3]*mee[4]*mee[4] - 2.0*mee[2]*mee[4]*mee[5] - mee[3]*mee[5]*mee[5] -
                2.0*mee[4]*mee[5]*cos(mee[6]) + onePaa*sin(mee[6]))
    dvzdk   = 2.0*sqrmup*invSs*invSs*(mee[3] + mee[3]*mee[4]*mee[4] - 2.0*mee[2]*mee[4]*mee[5] - mee[3]*mee[5]*mee[5] - 
                2.0*mee[4]*mee[5]*cos(mee[6]) + onePaa*sin(mee[6]))
            
    dvxdL   = -sqrmup*invSs*(onePaa*cos(mee[6]) + 2.0*mee[4]*mee[5]*sin(mee[6]))
    dvydL   = sqrmup*invSs*(-2.0*mee[4]*mee[5]*cos(mee[6]) - oneMaa*sin(mee[6]))
    dvzdL   = 2.0*sqrmup*invSs*(mee[5]*cos(mee[6]) - mee[4]*sin(mee[6]))  

    return SMatrix{6,6}(drxdp, drydp, drzdp, dvxdp, dvydp, dvzdp,
                        drxdf, drydf, drzdf, dvxdf, dvydf, dvzdf,
                        drxdg, drydg, drzdg, dvxdg, dvydg, dvzdg,
                        drxdh, drydh, drzdh, dvxdh, dvydh, dvzdh,
                        drxdk, drydk, drzdk, dvxdk, dvydk, dvzdk,
                        drxdL, drydL, drzdL, dvxdL, dvydL, dvzdL)
end