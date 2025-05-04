export eval_bihohmann
function eval_bihohmann(Δλ, tof, smaDrift; performDrift=true)
    Pg = 86164.0905
    mu = 3.986e5
    smaGeo = 42164.154046
    ωg = 360 / Pg

    smaHohmann = (smaDrift + smaGeo) / 2
    Ph = 2 * pi * sqrt(smaHohmann^3 / mu)
    Pd = 2 * pi * sqrt(smaDrift^3 / mu)

    lambdaDotHohmann = 360 / Ph - ωg
    lambdaDotDrift = 360 / Pd - ωg

    timeInDrift = tof - Ph
    if performDrift && timeInDrift < 0
        @warn "The Hohmann transfer sequence took too long, by $(-timeInDrift)"
        return 99999999999999999.9
    end
    Δλ_actual = lambdaDotHohmann * Ph
    if performDrift
        Δλ_actual += lambdaDotDrift * (timeInDrift)
    end
    return abs(Δλ - Δλ_actual)
end

export hohmann_cost
function hohmann_cost(r1::Real, r2::Real; mu=3.986e5)
    dv1 = sqrt(2 * mu * r2 / r1 / (r1 + r2)) - sqrt(mu / r1)
    dv2 = sqrt(mu / r2) - sqrt(2 * mu * r1 / r2 / (r1 + r2))
    return (abs(dv1) + abs(dv2), dv1, dv2)
end

export get_bihohmann_dv
function get_bihohmann_dv(Δλ, tof)
    # initial sma drift calculation, assuming whole time spent in drift orbit
    tofDays = tof / 86400.0
    rate = Δλ / tofDays
    smaOffsetDrift = abs(rate) * 78
    smaGeo = 42164.154046

    # Optimize around the actual drift sma needed
    optSol = Optim.optimize(x -> eval_bihohmann(Δλ, tof, x), smaGeo - 3 * smaOffsetDrift, smaGeo + 3 * smaOffsetDrift, GoldenSection())
    smaDrift = Optim.minimizer(optSol)

    (dv12, dv1, dv2) = hohmann_cost(smaGeo, smaDrift)
    (dv34, dv3, dv4) = hohmann_cost(smaDrift, smaGeo)

    return (dv12 + dv34)
end

export get_minimum_tof
function get_minimum_tof(Δλ)
    # initial sma drift calculation
    Pg = 86164.0905
    mu = 3.986e5
    smaGeo = 42164.154046
    ωg = 360 / Pg
    # Assume it takes about a day
    lambdaDot = Δλ / 86400.0
    ωh = lambdaDot + ωg
    Ph = 360.0 / ωh

    smaOffsetDrift = (Ph / 2 / pi)^2 * mu
    smaOffsetDrift = abs(smaOffsetDrift^(1 / 3) - smaGeo)

    # Optimize around the actual drift sma needed
    optSol = Optim.optimize(x -> eval_bihohmann(Δλ, 2 * pi * sqrt(x^3 / mu), x; performDrift=false), smaGeo - 3 * smaOffsetDrift, smaGeo + 3 * smaOffsetDrift, GoldenSection())
    smaDrift = Optim.minimizer(optSol)

    mintof = 2 * pi * sqrt(smaDrift^3 / mu)

    return mintof
end

export getnewvel
function getnewvel(x0, dv)
    pos = x0[1:3]
    vel = x0[4:6]
    r = unit(pos)
    c = unit(cross(r, vel))
    i = unit(cross(c, r))
    return vel + dv * i
end

export produce_bihohmann_results
function produce_bihohmann_results(Δλ)
    mu = 3.986e5
    minTof = get_minimum_tof(Δλ)
    sma_minTof = ((minTof / 2 / pi)^2 * mu)^(1 / 3)
    (dv12, dv1, dv2) = hohmann_cost(42164.0905, sma_minTof)

    times = Vector{Float64}()
    dvs = Vector{Float64}()

    push!(times, minTof)
    push!(dvs, 2 * abs(dv1))

    t = 86400.0 * ceil(minTof / 86400)

    while t <= 50 * 86400
        dv = get_bihohmann_dv(Δλ, t)

        push!(times, t)
        push!(dvs, dv)

        t += 86400.0
    end

    return (times, dvs)
end
