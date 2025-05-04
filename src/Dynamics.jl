function primerA(xIn)
    mu = 398600.4418
    rnorm = norm(xIn[1:3])
    rhat = 1 / rnorm * xIn[1:3]
    G = -mu / rnorm^3 * (I - 3 * rhat * rhat')

    # x = xIn[1]
    # y = xIn[2]
    # z = xIn[3]
    #
    # x2 = x * x
    # y2 = y * y
    # z2 = z * z
    #
    # denom = (x2 + y2 + z2)^(-5 / 2)
    #
    # G = zeros(3, 3)
    # G[1, 1] = mu * (2 * x2 - y2 - z2) * denom
    # G[1, 2] = 3 * mu * x * y * denom
    # G[1, 3] = 3 * mu * x * z * denom
    #
    # G[2, 1] = G[1, 2]
    # G[2, 2] = -mu * (x2 - 2 * y2 + z2) * denom
    # G[2, 3] = 3 * mu * y * z * denom
    #
    # G[3, 1] = G[1, 3]
    # G[3, 2] = G[2, 3]
    # G[3, 3] = -mu * (x2 + y2 - 2 * z2) * denom

    A = zeros(6, 6)
    A[1, 4] = 1.0
    A[2, 5] = 1.0
    A[3, 6] = 1.0
    A[4:6, 1:3] = G

    return A
end
function primerA_canonical(xIn)
    mu = 1.0
    rnorm = norm(xIn[1:3])
    rhat = 1 / rnorm * xIn[1:3]
    G = -mu / rnorm^3 * (I - 3 * rhat * rhat')

    A = zeros(6, 6)
    A[1, 4] = 1.0
    A[2, 5] = 1.0
    A[3, 6] = 1.0
    A[4:6, 1:3] = G

    return A
end

function dPrimer(x, p, t)
    A = primerA(x)
    return A * x
end

export dstate
function dstate(x, p, t)
    dx = zeros(6)
    r = norm(x[1:3])
    dx[1:3] = x[4:6]
    dx[4:6] = -398600.4418 / r^3 * x[1:3]

    return dx
end

function dstate_canonical(x, p, t)
    dx = zeros(6)
    r = norm(x[1:3])
    dx[1:3] = x[4:6]
    dx[4:6] = -1.0 / r^3 * x[1:3]

    return dx
end

function dxStm(x, p, t)
    state = x[1:6]
    stm = reshape(x[7:end], 6, 6)

    dx = zeros(size(x))
    r = norm(x[1:3])
    if r < 25
        dx[1:6] = dstate_canonical(state, p, t)
        A = primerA_canonical(state)
    else
        dx[1:6] = dstate(state, p, t)
        A = primerA(state)
    end
    dstm = A * stm
    dx[7:end] = reshape(dstm, 36, 1)

    return dx
end

function dxPrimer(x, p, t)
    state = x[1:6]
    primer = x[7:end]
    r = norm(x[1:3])
    if r < 25
        A = primerA_canonical(state)
        return [dstate_canonical(state, p, t); A * primer]
    else
        A = primerA(state)
        return [dstate(state, p, t); A * primer]
    end
end

export dstate_PS_ctrl
function dstate_PS_ctrl(x, p, t)
    ctrl = p[1]
    tof = p[2]
    tau = t * 2 / tof - 1
    N = size(ctrl, 1) - 1
    dx = zeros(6)
    r = norm(x[1:3])
    uaccel = testL2(tau, ctrl, N)
    dx[1:3] = x[4:6]
    dx[4:6] = -398600.4418 / r^3 * x[1:3] + uaccel

    return dx

end
function testPhi(t, l, N)
    tN, wN = _get_nodes_and_weights(N)
    phi = 1 / N / (N + 1) / Pl(tN[l], N) * (t^2 - 1) * dnPl(t, N, 1) / (t - tN[l])
    return phi
end
export testL2
function testL2(t, ctrl, N)
    lN = zeros(3)
    tN, _ = _get_nodes_and_weights(N)
    if t in tN
        ind = findall(in(tN), t)[1]
        return ctrl[ind, :]
    end

    for l in 1:N+1
        lN += ctrl[l, :] * testPhi(t, l, N)
    end
    return lN
end
function testPS()
    N = 20
    tN, _ = _get_nodes_and_weights(N)
    x0 = -1.0:0.01:1.0
    yN = sin.(2 * pi * tN)
    y0 = sin.(2 * pi * x0)
    y1 = [testL(t, x -> sin(2 * pi * x), N) for t in x0]

    set_theme!(theme_black())
    fig = Figure(size=(700, 600))

    ax = Axis(fig[1, 1])
    ax.xlabel = "t"
    ax.ylabel = "sin(t)"
    lines!(ax, x0, y0, color=:cyan, label="True")
    scatter!(ax, tN, yN, color=:orange, label="Knots")
    lines!(ax, x0, y1, color=:purple, label="PS")
    axislegend(ax, position=:cb)
    display(GLMakie.Screen(), fig)
end
