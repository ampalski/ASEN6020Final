export pseudospectral_impulsive_lambert
function pseudospectral_impulsive_lambert(
    x0::AbstractVector, #Starting state vector (position and velocity)
    x1::AbstractVector, #Final state vector (position or position and velocity)
    tof::AbstractFloat; #time of flight
    n::Int=-1, # Number of half revolutions between r1 and r2
    mu::AbstractFloat=3.986e5, #grav parameter for units used
    verbose::Bool=false, #whether or not to print debug statements
    poly_order::Int=20,
    constrain_u::Bool=false,
)
    POS_SCALING = 43000.0
    VEL_SCALING = 3.5
    # Input Checking
    if length(x0) != 6
        error("Initial state vector must be of length 6")
    end
    if length(x1) == 3
        positionOnly = true
    elseif length(x1) == 6
        positionOnly = false
    else
        error("Final state vector must be of length 3 or 6")
    end

    if n < 0
        n = _get_default_half_revs(0.5 * (norm(x0[1:3]) + norm(x1[1:3])), tof, mu=mu)
    end

    # DEFINE THE PROBLEM CONSTANTS
    r0 = x0[1:3]
    rf = x1[1:3]
    v0 = x0[4:6]
    vf = x1[4:6]
    verbose && (@info "Setting up initial seed")
    vxfer, vstop = basic_lambert(r0, v0, rf, tof, n, verbose=false, v2=vf)
    minaltsq = (DISTANCE_UNIT + 200)^2
    dv0 = vxfer - x0[4:6]
    dv1 = x1[4:6] - vstop
    dvinit = norm(dv0) + norm(dv1)
    maxDV = 2 * dvinit
    verbose && (@info "Initial Δv cost: $dvinit km/s")
    N = poly_order
    tN, wN = _get_nodes_and_weights(N)
    D = _differentiation_matrix(tN)
    nodes = Vector{PrimerNode}()
    push!(nodes, PrimerNode(0.0, x0[1:3], x0[4:6], vxfer))
    push!(nodes, PrimerNode(tof, x1[1:3], vstop, x1[4:6]))
    tauN = 0.5 .* (tof .* tN .+ tof)

    converged = false
    numBurns = 2


    while !converged && numBurns < 6

        # INITIALIZE THE MODEL
        verbose && (@info "Optimizing burn $numBurns")
        m = Model(Ipopt.Optimizer)
        if !verbose
            set_silent(m)
        end

        # INITIALIZE THE VARIABLES
        numSegments = length(nodes) - 1
        totalN = numSegments * (N + 1) + 2
        @variable(m, a1[1:totalN])
        @variable(m, a2[1:totalN])
        @variable(m, a3[1:totalN])
        @variable(m, a4[1:totalN])
        @variable(m, a5[1:totalN])
        @variable(m, a6[1:totalN])

        if constrain_u
            maxDV = 2 * dvinit
            @variable(m, -1.0 <= b1[1:numBurns] <= 1.0)
            @variable(m, -1.0 <= b2[1:numBurns] <= 1.0)
            @variable(m, -1.0 <= b3[1:numBurns] <= 1.0)
            # @variable(m, 0 <= b4[1:numBurns] <= tof)
            @variable(m, 0 <= b4[1:numBurns] <= 1.0)
        else
            @variable(m, b1[1:numBurns])
            @variable(m, b2[1:numBurns])
            @variable(m, b3[1:numBurns])
            # @variable(m, 0 <= b4[1:numBurns] <= tof)
            @variable(m, 0 <= b4[1:numBurns] <= 1.0)
        end

        # Constrain beginning and end

        fix(a1[1], r0[1] / POS_SCALING; force=true)
        fix(a2[1], r0[2] / POS_SCALING; force=true)
        fix(a3[1], r0[3] / POS_SCALING; force=true)
        fix(a4[1], v0[1] / VEL_SCALING; force=true)
        fix(a5[1], v0[2] / VEL_SCALING; force=true)
        fix(a6[1], v0[3] / VEL_SCALING; force=true)

        fix(a1[2], r0[1] / POS_SCALING; force=true)
        fix(a2[2], r0[2] / POS_SCALING; force=true)
        fix(a3[2], r0[3] / POS_SCALING; force=true)

        fix(a1[end], rf[1] / POS_SCALING; force=true)
        fix(a2[end], rf[2] / POS_SCALING; force=true)
        fix(a3[end], rf[3] / POS_SCALING; force=true)
        fix(a4[end], vf[1] / VEL_SCALING; force=true)
        fix(a5[end], vf[2] / VEL_SCALING; force=true)
        fix(a6[end], vf[3] / VEL_SCALING; force=true)

        fix(a1[end-1], rf[1] / POS_SCALING; force=true)
        fix(a2[end-1], rf[2] / POS_SCALING; force=true)
        fix(a3[end-1], rf[3] / POS_SCALING; force=true)

        fix(b4[1], 0.0; force=true)
        # fix(b4[end], tof; force=true)
        fix(b4[end], 1.0; force=true)
        @constraint(m, [i = 2:numBurns], b4[i] >= b4[i-1])

        # Set objective
        @objective(m, Min, sum(i -> sqrt(1e-14 + (b1[i]^2 + b2[i]^2 + b3[i]^2)), 1:numBurns))

        # Path constraint
        # @constraint(m, [k = 1:totalN], (a1[k] * POS_SCALING)^2 + (a2[k] * POS_SCALING)^2 + (a3[k] * POS_SCALING)^2 >= minaltsq)

        # Maneuver Starting Values
        for i in 1:numBurns
            dv = nodes[i].outVelocity - nodes[i].inVelocity
            set_start_value(b1[i], dv[1] / maxDV)
            set_start_value(b2[i], dv[2] / maxDV)
            set_start_value(b3[i], dv[3] / maxDV)
            if i != 1 && i != numBurns
                set_start_value(b4[i], nodes[i].time / tof)
            end
        end
        # Segment Endpoint Constraints
        @constraint(m, VEL_SCALING * a4[2] - VEL_SCALING * a4[1] - maxDV * b1[1] == 0.0)
        @constraint(m, VEL_SCALING * a5[2] - VEL_SCALING * a5[1] - maxDV * b2[1] == 0.0)
        @constraint(m, VEL_SCALING * a6[2] - VEL_SCALING * a6[1] - maxDV * b3[1] == 0.0)
        @constraint(m, VEL_SCALING * a4[end] - VEL_SCALING * a4[end-1] - maxDV * b1[end] == 0.0)
        @constraint(m, VEL_SCALING * a5[end] - VEL_SCALING * a5[end-1] - maxDV * b2[end] == 0.0)
        @constraint(m, VEL_SCALING * a6[end] - VEL_SCALING * a6[end-1] - maxDV * b3[end] == 0.0)
        if numSegments > 1
            for i in 2:numSegments
                ind = 1 + (i - 1) * (N + 1)
                @constraint(m, a1[ind+1] - a1[ind] == 0.0)
                @constraint(m, a2[ind+1] - a2[ind] == 0.0)
                @constraint(m, a3[ind+1] - a3[ind] == 0.0)
                @constraint(m, VEL_SCALING * a4[ind+1] - VEL_SCALING * a4[ind] - maxDV * b1[i] == 0.0)
                @constraint(m, VEL_SCALING * a5[ind+1] - VEL_SCALING * a5[ind] - maxDV * b2[i] == 0.0)
                @constraint(m, VEL_SCALING * a6[ind+1] - VEL_SCALING * a6[ind] - maxDV * b3[i] == 0.0)
            end
        end

        for i in 1:numSegments
            # PROVIDE STARTING VALUES
            x0p = [nodes[i].position...; nodes[i].outVelocity...]
            tau0 = nodes[i].time
            tauf = nodes[i+1].time
            prob = ODEProblem(dstate, x0p, (tau0, tauf), 0)
            sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
            guess_traj = [(t) -> sol(t)[i] for i in 1:6]
            tauN = 0.5 .* ((tauf - tau0) .* tN .+ (tauf + tau0))
            baseInd = (i - 1) * (N + 1) + 1

            for j in 1:N+1
                ind = baseInd + j
                tau = tauN[j]
                set_start_value(a1[ind], guess_traj[1](tau) / POS_SCALING)
                set_start_value(a2[ind], guess_traj[2](tau) / POS_SCALING)
                set_start_value(a3[ind], guess_traj[3](tau) / POS_SCALING)
                set_start_value(a4[ind], guess_traj[4](tau) / VEL_SCALING)
                set_start_value(a5[ind], guess_traj[5](tau) / VEL_SCALING)
                set_start_value(a6[ind], guess_traj[6](tau) / VEL_SCALING)
            end


            # Set dynamic constraints
            @constraint(m, [k = 1:N+1], (tof * b4[i+1] - tof * b4[i]) / 2 * VEL_SCALING * a4[baseInd+k] - sum(l -> D[k, l] * POS_SCALING * a1[baseInd+l], 1:N+1) == 0)

            @constraint(m, [k = 1:N+1], (tof * b4[i+1] - tof * b4[i]) / 2 * VEL_SCALING * a5[baseInd+k] - sum(l -> D[k, l] * POS_SCALING * a2[baseInd+l], 1:N+1) == 0)

            @constraint(m, [k = 1:N+1], (tof * b4[i+1] - tof * b4[i]) / 2 * VEL_SCALING * a6[baseInd+k] - sum(l -> D[k, l] * POS_SCALING * a3[baseInd+l], 1:N+1) == 0)

            @constraint(m, [k = 1:N+1], (tof * b4[i+1] - tof * b4[i]) / 2 * (-mu * POS_SCALING * a1[baseInd+k] / ((a1[baseInd+k] * POS_SCALING)^2 + (a2[baseInd+k] * POS_SCALING)^2 + (a3[baseInd+k] * POS_SCALING)^2)^(3 / 2)) - sum(l -> D[k, l] * VEL_SCALING * a4[baseInd+l], 1:N+1) == 0)

            @constraint(m, [k = 1:N+1], (tof * b4[i+1] - tof * b4[i]) / 2 * (-mu * POS_SCALING * a2[baseInd+k] / ((a1[baseInd+k] * POS_SCALING)^2 + (a2[baseInd+k] * POS_SCALING)^2 + (a3[baseInd+k] * POS_SCALING)^2)^(3 / 2)) - sum(l -> D[k, l] * VEL_SCALING * a5[baseInd+l], 1:N+1) == 0)

            @constraint(m, [k = 1:N+1], (tof * b4[i+1] - tof * b4[i]) / 2 * (-mu * POS_SCALING * a3[baseInd+k] / ((a1[baseInd+k] * POS_SCALING)^2 + (a2[baseInd+k] * POS_SCALING)^2 + (a3[baseInd+k] * POS_SCALING)^2)^(3 / 2)) - sum(l -> D[k, l] * VEL_SCALING * a6[baseInd+l], 1:N+1) == 0)
        end
        # SOLVE THE MODEL
        verbose && (@info "Solving...")
        optimize!(m)

        if termination_status(m) != LOCALLY_SOLVED &&
           termination_status(m) != ALMOST_LOCALLY_SOLVED
            @warn "Exiting early at $numBurns burns"
            return nodes
            error("Couldn't find solution.")
        end

        verbose && (@info "Looking to add a maneuver...")
        u1 = value.(b1) * maxDV
        u2 = value.(b2) * maxDV
        u3 = value.(b3) * maxDV
        u4 = value.(b4) * tof
        x_pos = value.(a1) * POS_SCALING
        y_pos = value.(a2) * POS_SCALING
        z_pos = value.(a3) * POS_SCALING
        vx = value.(a4) * VEL_SCALING
        vy = value.(a5) * VEL_SCALING
        vz = value.(a6) * VEL_SCALING
        for i in 1:numBurns
            index = (i - 1) * (N + 1) + 2
            if i > 1 && i < numBurns
                nodes[i].time = u4[i]
                nodes[i].position = [x_pos[index], y_pos[index], z_pos[index]]
            end
            if i > 1
                nodes[i].inVelocity = [vx[index-1], vy[index-1], vz[index-1]]
            end
            if i < numBurns
                nodes[i].outVelocity = [vx[index], vy[index], vz[index]]
            end
        end
        tfull, pmagfull = _get_full_primer_history(nodes)
        ind = argmax(pmagfull)
        pm = pmagfull[ind]
        if pm <= 1.01
            converged = true
            continue
        end

        numBurns += 1 # find the node it's in
        tmax = tfull[ind]
        baseInd = findlast([nodes[i].time for i in 1:length(nodes)] .< tmax)

        dv0 = nodes[baseInd].outVelocity - nodes[baseInd].inVelocity
        dv1 = nodes[baseInd+1].outVelocity - nodes[baseInd+1].inVelocity
        p0 = norm(dv0) == 0 ? zeros(3) : unit(dv0)
        p1 = norm(dv1) == 0 ? zeros(3) : unit(dv1)
        tofSegment = nodes[baseInd+1].time - nodes[baseInd].time
        x0p = [nodes[baseInd].position; nodes[baseInd].outVelocity]
        p0dot = _calculate_p0dot(p0, p1, tofSegment, x0p)
        stateHist = _get_full_state(p0, p0dot, x0p, tofSegment)
        tm, xm, drm, _ = _get_initial_mid_point(stateHist)

        # try here to just add the basic version and let ipopt take care of connecting

        # push!(nodes, PrimerNode(nodes[baseInd].time + tm, xm[1:3] + drm, xm[4:6], xm[4:6]))
        push!(nodes, PrimerNode(nodes[baseInd].time + tm, xm[1:3], xm[4:6], xm[4:6]))
        sort!(nodes)
        _correct_nodes!(nodes)
        # r_opt = [x_pos y_pos z_pos vx_pos vy_pos vz_pos]
        # return m, [u1 u2 u3], u4, r_opt
    end


    return nodes
end

function _correct_nodes!(nodes::Vector{PrimerNode})
    for i in 1:length(nodes)-1
        r0 = nodes[i].position
        v0 = nodes[i].inVelocity
        rf = nodes[i+1].position
        vf = nodes[i+1].outVelocity
        tof = nodes[i+1].time - nodes[i].time
        n = _get_default_half_revs(0.5 * (norm(r0) + norm(rf)), tof)
        vxfer, vstop = basic_lambert(r0, v0, rf, tof, n, verbose=false, v2=vf)
        nodes[i].outVelocity = vxfer
        nodes[i+1].inVelocity = vstop
    end
end
export produce_iPS_results
function produce_iPS_results(x0, x10)

    times = Vector{Float64}()
    dvs = Vector{Float64}()
    mu = 3.986e5

    t = 60 * 60.0

    while t <= 30 * 86400
        r, v = universalkepler(SA[x10[1:3]...], SA[x10[4:6]...], t, mu)
        x1 = [r; v]
        N = clamp(Int(floor(t / 3600)), 20, 100)
        tN, wN = _get_nodes_and_weights(N)
        nodes = pseudospectral_impulsive_lambert(x0, x1, t, poly_order=N, verbose=true, constrain_u=true, mu=mu)
        dv = 0.0
        for i in 1:length(nodes)
            dv += norm(nodes[i].outVelocity - nodes[i].inVelocity)
        end
        @info "t: $(t/86400), dv: $dv"

        push!(times, t)
        push!(dvs, dv)

        if t < 86400.0
            t += 6000.0
        elseif t < 5 * 86400.0
            t += 50000.0
        else
            t += 100000.0
        end
    end

    return (times, dvs)
end
