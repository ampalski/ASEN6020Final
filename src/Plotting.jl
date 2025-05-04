export plotFullPrimerHistory
function plotFullPrimerHistory(nodes::Vector{PrimerNode}; legend::Bool=false)
    tfull, pmagfull, impulses = _get_full_primer_history(nodes)
    set_theme!(theme_black())
    fig = Figure(size=(700, 600))

    ax = Axis(fig[1, 1])
    ax.xlabel = "Transfer Time (s)"
    ax.ylabel = "Primer Vector Magnitude"
    lines!(ax, tfull, pmagfull, color=:cyan, label="Primer History")
    scatter!(ax, impulses, ones(length(impulses)), color=:orange, label="Impulses")
    legend && axislegend(ax, position=:cb)
    display(GLMakie.Screen(), fig)
end

export plot_results
function plot_results(
    x0,
    x1,
    tof;
    primerNodes=nothing,
    collocationState=nothing,
    pseudoState=nothing,
)
    r0 = x0[1:3]
    rf = x1[1:3]
    v0 = x0[4:6]
    vf = x1[4:6]
    x1 = [x1...]
    mu = 3.986e5
    n = _get_default_half_revs(0.5 * (norm(x0[1:3]) + norm(x1[1:3])), tof, mu=mu)

    # initial orbit
    prob = ODEProblem(dstate, x0, (0, tof), 0)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    t = 0.0:10:tof
    initial = sol.(t)
    initialx = [initial[i][1] for i in 1:length(t)]
    initialy = [initial[i][2] for i in 1:length(t)]
    initialz = [initial[i][3] for i in 1:length(t)]

    # final orbit
    prob = ODEProblem(dstate, x1, (tof, 0), 0)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    final = sol.(t)
    finalx = [final[i][1] for i in 1:length(t)]
    finaly = [final[i][2] for i in 1:length(t)]
    finalz = [final[i][3] for i in 1:length(t)]

    vxfer, vstop = basic_lambert(r0, v0, rf, tof, n, verbose=false, v2=vf)
    x0p = [r0; vxfer]
    prob = ODEProblem(dstate, x0p, (0, tof), 0)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    lambert = sol.(t)
    lambertx = [lambert[i][1] for i in 1:length(t)]
    lamberty = [lambert[i][2] for i in 1:length(t)]
    lambertz = [lambert[i][3] for i in 1:length(t)]

    set_theme!(theme_black())
    fig = Figure(size=(800, 800))

    ax = Axis3(fig[1, 1])
    scatter!(ax, x0[1], x0[2], x0[3], color=:yellow, label="Initial Orbit")
    lines!(ax, initialx, initialy, initialz, color=:yellow)
    scatter!(ax, x1[1], x1[2], x1[3], color=:magenta, label="Final Orbit")
    lines!(ax, finalx, finaly, finalz, color=:magenta)

    lines!(ax, lambertx, lamberty, lambertz, color=:cyan, label="Lambert Soln")

    if !isnothing(primerNodes)
        for i in 1:length(primerNodes)-1
            x0p = [primerNodes[i].position...; primerNodes[i].outVelocity...]
            prob = ODEProblem(dstate, x0p, (primerNodes[i].time, primerNodes[i+1].time), 0)
            sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
            t = range(primerNodes[i].time, primerNodes[i+1].time, 100)
            lambert = sol.(t)
            lambertx = [lambert[i][1] for i in 1:length(t)]
            lamberty = [lambert[i][2] for i in 1:length(t)]
            lambertz = [lambert[i][3] for i in 1:length(t)]
            if i == 1
                lines!(ax, lambertx, lamberty, lambertz, color=:green, label="Primer")
            else
                lines!(ax, lambertx, lamberty, lambertz, color=:green)
            end
        end
    end

    if !isnothing(collocationState)
        lines!(ax, collocationState[:, 1], collocationState[:, 2], collocationState[:, 3], color=:orange, label="Collocation")
    end
    if !isnothing(pseudoState)
        lines!(ax, pseudoState[:, 1], pseudoState[:, 2], pseudoState[:, 3], color=:red, label="Pseudospectral")
    end
    axislegend(ax, position=:lb)

    xlims!(ax, -43000, 43000)
    ylims!(ax, -43000, 43000)
    # zlims!(ax, -43000, 43000)

    display(GLMakie.Screen(), fig)
end

export plot_optimality_chart
function plot_optimality_chart(results, Δλ)
    set_theme!(theme_black())
    fig = Figure(size=(700, 600))

    ax = Axis(fig[1, 1], yscale=log10)
    ax.xlabel = "Transfer Time (days)"
    ax.ylabel = "Δv (m/s)"
    # bihohmann
    if abs(Δλ) == 5
        t, dv = results["bihohmann5"]
    elseif abs(Δλ) == 10
        t, dv = results["bihohmann10"]
    elseif abs(Δλ) == 25
        t, dv = results["bihohmann25"]
    else
        error("Invalid change in longitude")
    end
    scatter!(ax, t ./ 86400, 1000 .* dv, label="Bi-Hohmann")

    # lambert
    if abs(Δλ) == 5
        t, dv = results["lambert5"]
    elseif abs(Δλ) == 10
        t, dv = results["lambert10"]
    elseif abs(Δλ) == 25
        t, dv = results["lambert25"]
    else
        error("Invalid change in longitude")
    end
    scatter!(ax, t ./ 86400, 1000 .* dv, label="Lambert")


    # primer
    if abs(Δλ) == 5
        t, dv = results["primer5"]
    elseif abs(Δλ) == 10
        t, dv = results["primer10"]
    elseif abs(Δλ) == 25
        t, dv = results["primer25"]
    else
        error("Invalid change in longitude")
    end
    scatter!(ax, t ./ 86400, 1000 .* dv, label="Primer Vectors")

    # PS
    if abs(Δλ) == 5
        t, dv = results["pseudo5"]
    elseif abs(Δλ) == 10
        t, dv = results["pseudo10"]
    elseif abs(Δλ) == 25
        t, dv = results["pseudo25"]
    else
        error("Invalid change in longitude")
    end
    scatter!(ax, t ./ 86400, 1000 .* dv, label="Pseudospectral")

    # iPS
    if abs(Δλ) == 5
        t, dv = results["ipseudo5"]
    elseif abs(Δλ) == 10
        t, dv = results["ipseudo10"]
    elseif abs(Δλ) == 25
        t, dv = results["ipseudo25"]
    else
        error("Invalid change in longitude")
    end
    scatter!(ax, t ./ 86400, 1000 .* dv, label="Impulsive Pseudospectral")

    axislegend(ax, position=:rt)
    display(GLMakie.Screen(), fig)
end

export plot_results2
function plot_results2(
    x0,
    x1,
    tof;
    primerNodes=nothing,
    pseudoState=nothing,
    control2=nothing,
    impPseudoNodes=nothing,
)
    r0 = x0[1:3]
    rf = x1[1:3]
    v0 = x0[4:6]
    vf = x1[4:6]
    x1 = [x1...]
    mu = 3.986e5
    n = _get_default_half_revs(0.5 * (norm(x0[1:3]) + norm(x1[1:3])), tof, mu=mu)

    # initial orbit
    prob = ODEProblem(dstate, x0, (0, tof), 0)
    solI = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    t = 0.0:10:tof
    # initial = sol.(t)
    # initialx = [initial[i][1] for i in 1:length(t)]
    # initialy = [initial[i][2] for i in 1:length(t)]
    # initialz = [initial[i][3] for i in 1:length(t)]
    # plotInitial = zeros(3, length(t))
    # for (ind, ti) in enumerate(t)
    #     plotInitial[:, ind] = ECI2ECF([initialx[ind], initialy[ind], initialz[ind]], ti)
    # end

    # final orbit
    prob = ODEProblem(dstate, x1, (tof, 0), 0)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    final = sol.(t)
    finalx = [final[i][1] for i in 1:length(t)]
    finaly = [final[i][2] for i in 1:length(t)]
    finalz = [final[i][3] for i in 1:length(t)]
    plotFinal = zeros(3, length(t))
    for (ind, ti) in enumerate(t)
        # plotFinal[:, ind] = ECI2ECF([finalx[ind], finaly[ind], finalz[ind]], ti)
        plotFinal[:, ind] = eci2ric(solI(ti), [finalx[ind], finaly[ind], finalz[ind]])
    end

    vxfer, vstop = basic_lambert(r0, v0, rf, tof, n, verbose=false, v2=vf)
    x0p = [r0; vxfer]
    prob = ODEProblem(dstate, x0p, (0, tof), 0)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    lambert = sol.(t)
    lambertx = [lambert[i][1] for i in 1:length(t)]
    lamberty = [lambert[i][2] for i in 1:length(t)]
    lambertz = [lambert[i][3] for i in 1:length(t)]
    plotlambert = zeros(3, length(t))
    for (ind, ti) in enumerate(t)
        # plotlambert[:, ind] = ECI2ECF([lambertx[ind], lamberty[ind], lambertz[ind]], ti)
        plotlambert[:, ind] = eci2ric(solI(ti), [lambertx[ind], lamberty[ind], lambertz[ind]])
    end

    x10 = [finalx[1], finaly[1], finalz[1]]
    Δλ = anglevec(x0[1:3], x10) * 180 / pi
    temp = cross(x0[1:3], x10)
    if temp[3] < 0
        Δλ *= -1
    end
    tofDays = tof / 86400.0
    rate = Δλ / tofDays
    smaOffsetDrift = abs(rate) * 78
    smaGeo = 42164.154046
    # Optimize around the actual drift sma needed
    optSol = Optim.optimize(x -> eval_bihohmann(Δλ, tof, x), smaGeo - 3 * smaOffsetDrift, smaGeo + 3 * smaOffsetDrift, GoldenSection())
    smaDrift = Optim.minimizer(optSol)
    driftsma = (smaDrift + smaGeo) / 2
    Pxfer = pi * sqrt(driftsma^3 / mu)
    Pdrift = tof - 2 * Pxfer
    (dv12, dv1, dv2) = hohmann_cost(smaGeo, smaDrift)
    (dv34, dv3, dv4) = hohmann_cost(smaDrift, smaGeo)
    v0p = getnewvel(x0, dv1)
    prob = ODEProblem(dstate, vcat(r0, v0p), (0, Pxfer), 0)
    sol1 = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    rv1 = sol1[end]
    v1p = getnewvel(rv1, dv2)
    prob = ODEProblem(dstate, vcat(rv1[1:3], v1p), (Pxfer, Pxfer + Pdrift), 0)
    sol2 = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    rv2 = sol2[end]
    v2p = getnewvel(rv2, dv3)
    prob = ODEProblem(dstate, vcat(rv2[1:3], v2p), (Pxfer + Pdrift, tof), 0)
    sol3 = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    plotbihohmann = zeros(3, length(t))
    for (ind, ti) in enumerate(t)
        if ti <= Pxfer
            temp = sol1(ti)
        elseif ti <= (Pxfer + Pdrift)
            temp = sol2(ti)
        else
            temp = sol3(ti)
        end
        # plotbihohmann[:, ind] = ECI2ECF(temp[1:3], ti)
        plotbihohmann[:, ind] = eci2ric(solI(ti), temp[1:3])
    end


    set_theme!(theme_black())
    fig = Figure(size=(800, 800))

    ax = Axis3(fig[1, 1])
    # x0plot = ECI2ECF(x0[1:3], 0.0)
    # x1plot = ECI2ECF(x1[1:3], tof)
    x0plot = zeros(3)
    x1plot = eci2ric(solI[end], x1[1:3])
    scatter!(ax, x0plot[1], x0plot[2], x0plot[3], color=:yellow, label="Initial Orbit")
    # lines!(ax, plotInitial[1, :], plotInitial[2, :], plotInitial[3, :], color=:yellow)
    scatter!(ax, x1plot[1], x1plot[2], x1plot[3], color=:magenta, label="Final Orbit")
    lines!(ax, plotFinal[1, :], plotFinal[2, :], plotFinal[3, :], color=:magenta)

    lines!(ax, plotlambert[1, :], plotlambert[2, :], plotlambert[3, :], color=:cyan, label="Lambert Soln")
    lines!(ax, plotbihohmann[1, :], plotbihohmann[2, :], plotbihohmann[3, :], label="BiHohmann Transfer")

    if !isnothing(primerNodes)
        for i in 1:length(primerNodes)-1
            x0p = [primerNodes[i].position...; primerNodes[i].outVelocity...]
            prob = ODEProblem(dstate, x0p, (primerNodes[i].time, primerNodes[i+1].time), 0)
            sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
            t = range(primerNodes[i].time, primerNodes[i+1].time, 100)
            lambert = sol.(t)
            lambertx = [lambert[i][1] for i in 1:length(t)]
            lamberty = [lambert[i][2] for i in 1:length(t)]
            lambertz = [lambert[i][3] for i in 1:length(t)]
            plotlambert = zeros(3, length(t))
            for (ind, ti) in enumerate(t)
                # plotlambert[:, ind] = ECI2ECF([lambertx[ind], lamberty[ind], lambertz[ind]], ti)
                plotlambert[:, ind] = eci2ric(solI(ti), [lambertx[ind], lamberty[ind], lambertz[ind]])
            end
            if i == 1
                lines!(ax, plotlambert[1, :], plotlambert[2, :], plotlambert[3, :], color=:green, label="Primer")
            else
                lines!(ax, plotlambert[1, :], plotlambert[2, :], plotlambert[3, :], color=:green)
            end
        end
    end

    if !isnothing(pseudoState)
        N = size(pseudoState, 1) - 1
        plotPseudo = zeros(3, N + 1)
        tN, wN = _get_nodes_and_weights(N)
        tauN = (tof .* tN .+ tof) ./ 2
        for (ind, ti) in enumerate(tauN)
            # plotPseudo[:, ind] = ECI2ECF([pseudoState[ind, 1], pseudoState[ind, 2], pseudoState[ind, 3]], ti)
            plotPseudo[:, ind] = eci2ric(solI(ti), [pseudoState[ind, 1], pseudoState[ind, 2], pseudoState[ind, 3]])
        end
        lines!(ax, plotPseudo[1, :], plotPseudo[2, :], plotPseudo[3, :], color=:red, label="Pseudospectral")
    end
    if !isnothing(control2)
        prob = ODEProblem(dstate_PS_ctrl, x0, (0, tof), (control2, tof))
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        final = sol.(t)
        finalx = [final[i][1] for i in 1:length(t)]
        finaly = [final[i][2] for i in 1:length(t)]
        finalz = [final[i][3] for i in 1:length(t)]
        plotFinal = zeros(3, length(t))
        for (ind, ti) in enumerate(t)
            # plotFinal[:, ind] = ECI2ECF([finalx[ind], finaly[ind], finalz[ind]], ti)
            plotFinal[:, ind] = eci2ric(solI(ti), [finalx[ind], finaly[ind], finalz[ind]])
        end
        lines!(ax, plotFinal[1, :], plotFinal[2, :], plotFinal[3, :], color=:green, label="Propagated PS Control")
    end
    if !isnothing(impPseudoNodes)
        for i in 1:length(impPseudoNodes)-1
            x0p = [impPseudoNodes[i].position...; impPseudoNodes[i].outVelocity...]
            prob = ODEProblem(dstate, x0p, (impPseudoNodes[i].time, impPseudoNodes[i+1].time), 0)
            sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
            t = range(impPseudoNodes[i].time, impPseudoNodes[i+1].time, 100)
            lambert = sol.(t)
            lambertx = [lambert[i][1] for i in 1:length(t)]
            lamberty = [lambert[i][2] for i in 1:length(t)]
            lambertz = [lambert[i][3] for i in 1:length(t)]
            plotlambert = zeros(3, length(t))
            for (ind, ti) in enumerate(t)
                # plotlambert[:, ind] = ECI2ECF([lambertx[ind], lamberty[ind], lambertz[ind]], ti)
                plotlambert[:, ind] = eci2ric(solI(ti), [lambertx[ind], lamberty[ind], lambertz[ind]])
            end
            if i == 1
                lines!(ax, plotlambert[1, :], plotlambert[2, :], plotlambert[3, :], color=:salmon, label="Impulsive Pseudospectral")
            else
                lines!(ax, plotlambert[1, :], plotlambert[2, :], plotlambert[3, :], color=:salmon)
            end
        end
    end
    axislegend(ax, position=:lb)
    ax.xlabel = "Radial"
    ax.ylabel = "In-Track"
    ax.zlabel = "Cross-Track"

    # xlims!(ax, -43000, 43000)
    # ylims!(ax, -43000, 43000)
    # zlims!(ax, -43000, 43000)

    display(GLMakie.Screen(), fig)
end

export ECI2ECF
function ECI2ECF(pos, t)
    gmst = GMST(t)
    ecf2eci = zeros(3, 3)
    ecf2eci[3, 3] = 1.0
    c = cos(gmst)
    s = sin(gmst)
    ecf2eci[1, 1] = c
    ecf2eci[2, 2] = c
    ecf2eci[1, 2] = -s
    ecf2eci[2, 1] = s

    return ecf2eci' * pos
end
export GMST
function GMST(t)
    baseEpoch = 2460243.077199
    epoch = baseEpoch + t / 86400.0

    JD0 = floor(epoch) + 0.5
    JD0 > epoch && (JD0 -= 1)
    H = 24.0 * (epoch - JD0)

    D = epoch - 2451545.0
    T = D / 36525
    GMST = mod(6.697375 + 0.065709824279 * D + 1.0027379 * H + 0.0854103 * T + 0.0000258 * T^2, 24.0)
    GMST = GMST * pi / 12 #to degrees, to radians
    return GMST
end
export eci2ric
function eci2ric(target::AbstractArray, chase::AbstractArray)
    R = unit(target[1:3])
    C = unit(cross(R, target[4:6]))
    I = unit(cross(C, R))

    posdiff = chase[1:3] - target[1:3]

    ricpos = [R'; I'; C'] * posdiff
    tgtR = norm(target[1:3])
    chaseR = norm(chase[1:3])
    sph = zeros(3)

    sph[1] = chaseR - tgtR
    sph[2] = tgtR * atan(ricpos[2], tgtR + ricpos[1])
    sph[3] = tgtR * asin(ricpos[3] / chaseR)
    return sph
end
