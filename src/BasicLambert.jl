include("LambertBattin.jl")
include("LambertGooding.jl")
include("LambertIzzo.jl")
include("ParticleSwarm.jl")

export lambert_battin, lambert_gooding

#Tries several attempts with n and pick to find the min delta-v case
function basic_lambert(
    r1::AbstractVector, #Starting position vector
    v1::AbstractVector, #Starting velocity vector
    r2::AbstractVector, #Ending position vector
    tof::AbstractFloat, #time of flight
    n::Int; # Number of half revolutions between r1 and r2
    mu::AbstractFloat=3.986e5, #grav parameter for units used
    verbose::Bool=false, #whether or not to print debug statements
    v2::AbstractVector=zeros(3),
    runPSO::Bool=false,
)
    includestop = norm(v2) == 0 ? false : true
    l = 18
    vx, vs = lambert_izzo(r1, r2, tof, mu)
    l = 18 + length(vx)
    vxfer = zeros(3, l)
    vstop = zeros(3, l)

    dvs = zeros(l)
    nrev = Int(floor(n / 2))
    vxfer[:, 1], vstop[:, 1] = lambert_gooding(r1, r2, tof, n, pick=1, mu=mu)
    vxfer[:, 2], vstop[:, 2] = lambert_gooding(r1, r2, tof, n, pick=2, mu=mu)
    vxfer[:, 3], vstop[:, 3] = lambert_gooding(r1, r2, tof, max(n - 1, 0), pick=1, mu=mu)
    vxfer[:, 4], vstop[:, 4] = lambert_gooding(r1, r2, tof, max(n - 1, 0), pick=2, mu=mu)
    vxfer[:, 5], vstop[:, 5] = lambert_gooding(r1, r2, tof, n + 1, pick=1, mu=mu)
    vxfer[:, 6], vstop[:, 6] = lambert_gooding(r1, r2, tof, n + 1, pick=2, mu=mu)
    vxfer[:, 7], vstop[:, 7] = lambert_battin(r1, v1, r2, tof, nrev, 1, 1, mu=mu)
    vxfer[:, 8], vstop[:, 8] = lambert_battin(r1, v1, r2, tof, nrev, 1, -1, mu=mu)
    vxfer[:, 9], vstop[:, 9] = lambert_battin(r1, v1, r2, tof, nrev, -1, 1, mu=mu)
    vxfer[:, 10], vstop[:, 10] = lambert_battin(r1, v1, r2, tof, nrev, -1, -1, mu=mu)
    vxfer[:, 11], vstop[:, 11] = lambert_battin(r1, v1, r2, tof, max(nrev - 1, 0), 1, 1, mu=mu)
    vxfer[:, 12], vstop[:, 12] = lambert_battin(r1, v1, r2, tof, max(nrev - 1, 0), 1, -1, mu=mu)
    vxfer[:, 13], vstop[:, 13] = lambert_battin(r1, v1, r2, tof, max(nrev - 1, 0), -1, 1, mu=mu)
    vxfer[:, 14], vstop[:, 14] = lambert_battin(r1, v1, r2, tof, max(nrev - 1, 0), -1, -1, mu=mu)
    vxfer[:, 15], vstop[:, 15] = lambert_battin(r1, v1, r2, tof, nrev + 1, 1, 1, mu=mu)
    vxfer[:, 16], vstop[:, 16] = lambert_battin(r1, v1, r2, tof, nrev + 1, 1, -1, mu=mu)
    vxfer[:, 17], vstop[:, 17] = lambert_battin(r1, v1, r2, tof, nrev + 1, -1, 1, mu=mu)
    vxfer[:, 18], vstop[:, 18] = lambert_battin(r1, v1, r2, tof, nrev + 1, -1, -1, mu=mu)

    ind = 19
    for i in eachindex(vs)
        vxfer[:, ind] = vx[i]
        vstop[:, ind] = vs[i]
        ind += 1
    end

    for i in 1:l
        if includestop
            dv = norm(vxfer[:, i] - v1) + norm(v2 - vstop[:, i])
        else
            dv = norm(vxfer[:, i] - v1)
        end
        dvs[i] = vxfer[:, i] == zeros(3) ? 9999999.0 : dv

        rf, vf = universalkepler(SA[r1...], SA[vxfer[:, i]...], tof, mu)

        if norm(rf - r2) > 25 && norm(rf - r2) <= 50

            dvs[i] = 9999999999999.0
        elseif norm(rf - r2) > 50
            dvs[i] = 9999999999999.0
            # verbose && println("i is $i, missed by $(norm(rf-r2))")
        end
    end
    dvs[isnan.(dvs)] .= 99999999999999999.0
    dvmin, ind = findmin(dvs)
    verbose && @info "Lambert solver chose i=$ind with dv=$dvmin"
    if dvmin > 50 && runPSO
        verbose && println("Attempting to solve via PSO...")
        v1PSO = ParticleSwarmOptimizer([[v1[1] - 1, v1[1] + 1], [v1[2] - 1, v1[2] + 1], [v1[3] - 1, v1[3] + 1]], x -> LambertPSOEval(r1, r2, tof, x), verbose=false)
        rf, vf = universalkepler(SA[r1...], SA[v1PSO...], tof, mu)
        if norm(rf - r2) > 10
            verbose && println("Opening up PSO...")
            v1PSO = ParticleSwarmOptimizer([[v1[1] - 2, v1[1] + 2], [v1[2] - 2, v1[2] + 2], [v1[3] - 2, v1[3] + 2]], x -> LambertPSOEval(r1, r2, tof, x), verbose=false, maxIter=1000)
            rf, vf = universalkepler(SA[r1...], SA[v1PSO...], tof, mu)
        end
        norm(rf - r2) > 10 && error("No valid solution found\n$r1\n$r2\n$v1\n$tof")
        verbose && println("PSO Solution with dv=$(norm(v1PSO-v1))")
        return (v1PSO, vf)
    end


    return (vxfer[:, ind], vstop[:, ind])
end

function LambertPSOEval(x0, xf, dt, vxfer)
    rf, _ = universalkepler(SA[x0...], SA[vxfer...], dt, 3.986e5)

    return norm(xf - rf)
end

export produce_lambert_results
function produce_lambert_results(x0, x10)

    times = Vector{Float64}()
    dvs = Vector{Float64}()
    mu = 3.986e5

    t = 60 * 60.0

    while t <= 50 * 86400
        r, v = universalkepler(SA[x10[1:3]...], SA[x10[4:6]...], t, mu)
        x1 = [r; v]
        n = _get_default_half_revs(0.5 * (norm(x0[1:3]) + norm(x1[1:3])), t, mu=mu)
        vxfer, vstop = basic_lambert(x0[1:3], x0[4:6], x1[1:3], t, n, mu=mu, verbose=false, v2=x1[4:6])
        dv = norm(vxfer - x0[4:6]) + norm(x1[4:6] - vstop)

        push!(times, t)
        push!(dvs, dv)

        if t < 86400.0
            t += 4000.0
        elseif t < 5 * 86400.0
            t += 40000.0
        else
            t += 80000.0
        end
    end

    return (times, dvs)
end
