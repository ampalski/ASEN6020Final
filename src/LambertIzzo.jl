
export lambert_izzo
function lambert_izzo(
    r1::AbstractVector, #Starting position vector
    r2::AbstractVector, #Ending position vector
    tof::AbstractFloat, #time of flight
    mu::AbstractFloat=3.986e5, #grav parameter for units used
)
    if tof < 0 || mu < 0
        error("Invalid time of flight or gravitational parameter")
    end

    c = r2 - r1
    r1n = norm(r1)
    r2n = norm(r2)
    cn = norm(c)
    s = 0.5 * (r1n + r2n + cn)

    ir1 = unit(r1)
    ir2 = unit(r2)
    ih = cross(ir1, ir2)
    ih = unit(ih)

    λ2 = 1 - cn / s
    λ = sqrt(λ2)

    if (r1[1] * r2[2] - r1[2] * r2[1]) < 0
        λ = -λ
        it1 = cross(ir1, ih)
        it2 = cross(ir2, ih)
    else
        it1 = cross(ih, ir1)
        it2 = cross(ih, ir2)
    end

    T = sqrt(2 * mu / s^3) * tof
    gamma = sqrt(mu * s / 2)
    rho = (r1n - r2n) / cn
    sigma = sqrt(1 - rho^2)

    x_all = _findxy(λ, T)

    v1 = Vector{Vector{Float64}}()
    v2 = Vector{Vector{Float64}}()
    for x in x_all
        y = sqrt(1.0 - λ2 + λ2 * x * x)
        vr1 = gamma * ((λ * y - x) - rho * (λ * y + x)) / r1n
        vr2 = -gamma * ((λ * y - x) + rho * (λ * y + x)) / r2n
        vt = gamma * sigma * (y + λ * x)
        vt1 = vt / r1n
        vt2 = vt / r2n
        push!(v1, (vr1 * ir1 + vt1 * it1))
        push!(v2, (vr2 * ir2 + vt2 * it2))
    end
    return (v1, v2)
end

function _findxy(λ, T)
    λ2 = λ * λ
    λ3 = λ2 * λ

    Nmax = floor(T / pi)

    T00 = acos(λ) + λ * sqrt(1 - λ2)
    T0 = T00 + Nmax * pi
    T1 = 2 / 3 * (1 - λ3)
    if T < T0 && Nmax > 0
        iter = 0
        errval = 1.0
        Tmin = T0
        x_old = 0.0
        x_new = 0.0
        while true
            dT, d2T, d3T = dTdx(x_old, Tmin, λ)
            x_new = x_old - dT * d2T / (d2T * d2T - dT * d3T / 2)
            errval = abs(x_old - x_new)
            if errval < 1e-13 || iter > 13
                break
            end
            Tmin = x2tof(x_new, Nmax, λ)
            x_old = x_new
            iter += 1
        end
        if Tmin > T
            Nmax -= 1
        end
    end

    # Start the solution vectors, Find initial x0
    x_all = Vector{Float64}()
    x_cand = 0.0

    if T >= T0
        x_cand = -(T - T00) / (T - T00 + 4)
    elseif T < T1
        x_cand = T1 * (T1 - T) / (2.0 / 5.0 * (1 - λ2 * λ3) * T) + 1
    else
        x_cand = (T / T00)^(0.69314718055994529 / log(T1 / T00)) - 1
    end

    # Iterate
    x_cand = householder(T, x_cand, 0, 1e-5, 15, λ)
    push!(x_all, x_cand)

    # Search multi-rev
    for i in 1:Nmax
        temp = ((i * pi + pi) / (8 * T))^(2.0 / 3)
        xl = (temp - 1) / (temp + 1)
        xl = householder(T, xl, i, 1e-8, 15, λ)

        temp = ((8 * T) / (i * pi))^(2.0 / 3)
        xr = (temp - 1) / (temp + 1)
        xr = householder(T, xr, i, 1e-8, 15, λ)

        push!(x_all, xl)
        push!(x_all, xr)
    end

    return x_all
end

export householder
function householder(T, x, N, eps, iter_max, λ)
    x0 = x
    iter = 0
    errval = 1.0
    xnew = 0.0
    tof = 0.0
    delta = 0.0
    dT = 0.0
    d2T = 0.0
    d3T = 0.0
    while errval > eps && iter < iter_max
        tof = x2tof(x0, N, λ)
        dT, d2T, d3T = dTdx(x0, tof, λ)
        delta = tof - T
        dT2 = dT * dT
        xnew = x0 - delta * (dT2 - delta * d2T / 2.0) / (dT * (dT2 - delta * d2T) + d3T * delta * delta / 6.0)
        errval = abs(x0 - xnew)
        x0 = xnew
        iter += 1
    end
    return xnew
end

export x2tof, x2tof2
function x2tof2(x, N, λ)
    a = 1 / (1 - x * x)

    if a > 0
        alpha = 2 * acos(x)
        beta = 2 * asin(sqrt(λ * λ / a))
        if λ < 0
            beta = -beta
        end
        T = ((a * sqrt(a) * ((alpha - sin(alpha)) - (beta - sin(beta)) + 2 * pi * N)) / 2)
        return T
    else
        alpha = 2 * acosh(x)
        beta = 2 * asinh(-λ * λ / a)
        if λ < 0
            beta = -beta
        end
        T = (-a * sqrt(-a) * ((beta - sinh(beta)) - (alpha - sinh(alpha))) / 2)
        return T
    end
end

function x2tof(x, N, λ)
    battin = 0.01
    lagrange = 0.2
    dist = abs(x - 1)
    if dist < lagrange && dist > battin
        return x2tof2(x, N, λ)
    end

    K = λ * λ
    E = x * x - 1
    rho = abs(E)
    z = sqrt(1 + K * E)

    if dist < battin
        eta = z - λ * x
        S1 = 0.5 * (1 - λ - x * eta)
        Q = hypergeometricF(S1, 1e-11)
        T = (eta * eta * eta * Q + 4 * λ * eta) / 2.0 + N * pi / rho^(1.5)
        return T
    else
        y = sqrt(rho)
        g = x * z - λ * E
        d = 0.0
        if E < 0
            l = acos(g)
            d = N * pi + l
        else
            f = y * (z - λ * x)
            d = log(f + g)
        end
        T = (x - λ * z - d / y) / E
        return T
    end
end

function dTdx(x, T, λ)
    λ2 = λ * λ
    λ3 = λ2 * λ
    umx2 = 1 - x^2
    umx2i = 1 / umx2
    y = sqrt(1 - λ2 * umx2)
    y2 = y * y
    y3 = y2 * y

    dT = umx2i * (3 * T * x - 2 + 2 * λ3 * x / y)
    d2T = umx2i * (3 * T + 5 * x * dT + 2 * (1 - λ2) * λ3 / y3)
    d3T = umx2i * (7 * x * d2T + 8 * dT - 6 * (1 - λ2) * λ2 * λ3 * x / y3 / y2)

    return (dT, d2T, d3T)
end

function hypergeometricF(z, tol)
    Sj = 1.0
    Cj = 1.0
    errval = 1.0
    Cj1 = 0.0
    Sj1 = 0.0
    j = 0
    while errval > tol
        Cj1 = Cj * (3 + j) * (1 + j) / (2.5 + j) * z / (j + 1)
        Sj1 = Sj + Cj1
        errval = abs(Cj1)
        Sj = Sj1
        Cj = Cj1
        j = j + 1
    end
    return Sj
end
