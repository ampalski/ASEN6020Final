export universalkepler
"""
    (newPosition, newVelocity) = universalkepler(r0, v0, dt, mu)

Propagate an orbit state according to 2-body dynamics

# Arguments
- `r0`: The initial position vector (3x1 vector)
- `v0`: The initial velocity vector (3x1 vector)
- `dt`: Amount of time the state should be propagated over
- `mu`: Gravitational parameter, with units consistent with the other inputs
"""
function universalkepler(
    r0::SVector{3,Float64},
    v0::SVector{3,Float64},
    dt::Float64,
    mu::Float64,
)

    r0n = norm(r0)
    v0n = norm(v0)
    temp = dot(r0, v0) / sqrt(mu)

    alpha = -v0n * v0n / mu + 2 / r0n

    if alpha > 0.000001 && abs(alpha - 1) > 0.000001 # circle or ellipse
        ChiOld = sqrt(mu) * dt * alpha
    elseif alpha > 0.000001 # circle or ellipse with bad convergence
        ChiOld = sqrt(mu) * dt * alpha * 0.97
    elseif abs(alpha) < 0.000001 #parabolic
        h = cross(r0, v0)
        hn = norm(h)
        p = hn * hn / mu
        s = 0.5 * acot(3 * sqrt(mu / (p^3)) * dt)
        w = atan((tan(s))^(1 / 3))
        ChiOld = sqrt(p) * 2 * cot(2 * w)
        alpha = 0
    else #hyperbolic
        a = 1 / alpha
        ChiOld = sign(dt) * sqrt(-a)
        ChiOld *= log(-2 * mu * alpha * dt /
                      (dot(r0, v0) + sign(dt) * sqrt(-mu * a) * (1 - r0n * alpha)))
    end

    Chi = 999.0
    ktr = 0
    c2 = 0.0
    c3 = 0.0
    Psi = 0.0
    r = 0.0

    while true
        Psi = ChiOld * ChiOld * alpha
        (c2, c3) = _findc2c3(Psi)
        r = ChiOld^2 * c2 + temp * ChiOld * (1 - Psi * c3) + r0n * (1 - Psi * c2)
        dChir = sqrt(mu) * dt - ChiOld^3 * c3 -
                temp * ChiOld^2 * c2 - r0n * ChiOld * (1 - Psi * c3)
        Chi = ChiOld + dChir / r

        if abs(Chi - ChiOld) < 1e-6 || ktr > 100
            break
        end

        ChiOld = Chi
        ktr += 1
    end

    Chisq = Chi * Chi
    f = 1 - Chisq / r0n * c2
    g = dt - Chisq * Chi * c3 / sqrt(mu)
    gd = 1 - Chisq / r * c2
    fd = sqrt(mu) / r / r0n * Chi * (Psi * c3 - 1)

    r = f * r0 + g * v0
    v = fd * r0 + gd * v0

    return (r, v)

end

function _findc2c3(Psi)

    if Psi > 1e-6
        rtPsi = sqrt(Psi)
        c2 = (1 - cos(rtPsi)) / Psi
        c3 = (rtPsi - sin(rtPsi)) / sqrt(Psi^3)
    else
        if Psi < -1e-6
            rtPsi = sqrt(-Psi)
            c2 = (1 - cosh(rtPsi)) / Psi
            c3 = (sinh(rtPsi) - rtPsi) / sqrt((-Psi)^3)
        else
            c2 = 0.5
            c3 = 1 / 6
        end
    end

    return (c2, c3)
end

