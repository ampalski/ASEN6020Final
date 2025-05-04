function a1(eta, n)
    return eta * n * n / ((2 * n)^2 - 1)
end

function a2(U, n)
    if isodd(n)
        return U * 2 * (3 * n + 1) * (6 * n - 1) / 9 / (4 * n - 1) / (4 * n + 1)
    else
        return U * 2 * (3 * n + 2) * (6 * n + 1) / 9 / (4 * n + 3) / (4 * n + 1)
    end
end

function b1(n)
    return 1.0
end


function _LambertMinTOF(r1, r2, n, dir, mu)
    r1n = norm(r1)
    r2n = norm(r2)
    c = norm(r2 - r1)
    s = (r1n + r2n + c) / 2
    an = 1.001 * s / 2

    fa = 20.0
    ξ = 20.0
    η = 20.0
    itr = 0
    while abs(fa) < 1e-16
        alpha = 2 * asin(sqrt(s / 2 / an))
        beta = 2 * dir * asin(sqrt((s - c) / 2 / an))

        ξ = alpha - beta
        η = sin(alpha) - sin(beta)

        fa = ((6 * n * pi + 3 * ξ - η) * (sin(ξ) + η) - 8 * (1 - cos(ξ)))
        fpa = ((6 * n * pi + 3 * ξ - η) * (cos(ξ) + cos(alpha)) +
               (3 - cos(alpha)) * (sin(ξ) + η) - 8 * sin(ξ)) * (-1 / an * tan(alpha / 2)) +
              ((6 * n * pi + 3 * ξ - η) * (-cos(ξ) - cos(beta)) +
               (-3.0 + cos(beta)) * (sin(ξ) + η) + 8 * sin(ξ)) * (-1 / an * tan(beta / 2))

        an = an - fa / fpa
        itr += 1
        if itr > 1000
            error("too many iterations")
        end
    end
    return sqrt(an^3 / mu) * (2 * n * pi + ξ - η)

end

# Following Vallado Version 5 Alg 61
# Consider adding Alg 57 in the Lambert fn to check for tofs/nrevs that make sense
# And 58 to make sure they don't hit the earth
function lambert_battin(
    r1::AbstractVector, #Starting position vector
    v1::AbstractVector, #Starting velocity vector
    r2::AbstractVector, #Ending position vector
    tof::AbstractFloat, #time of flight
    n::Int, # Number of revolutions between r1 and r2
    dir::Int, # Direction of motion (+1 for short, -1 for long)
    energy::Int; # +1 for low energy case, -1 for high energy case
    mu::AbstractFloat=3.986e5 #grav parameter for units used
)

    if _LambertMinTOF(r1, r2, n, dir, mu) > tof
        return (zeros(3), zeros(3))
    end

    # Set up initial geometric values
    r1n = norm(r1)
    r2n = norm(r2)

    cdθ = dot(r1, r2) / r1n / r2n
    sdθ = dir * sqrt(1 - cdθ^2)
    dθ = atan(sdθ, cdθ)
    if dθ < 0
        dθ += 2.0 * pi
    end

    c = norm(r2 - r1)
    s = (r1n + r2n + c) / 2

    λ = 1 / s * sqrt(r1n * r2n) * cos(dθ / 2)
    l = ((1 - λ) / (1 + λ))^2
    m = 8 * mu * tof^2 / s^3 / (1 + λ)^6

    a = 0.0
    p = 0.0
    e = 0.0

    # High energy case
    if energy < 0 && n > 0
        x = 1e-20
        xold = 10.0
        loops = 0
        while abs(x - xold) > 1e-16 && loops < 30
            sqrtx = sqrt(x)
            h1 = (l + x) * (1 + 2 * x + l) / 2 / (l - x^2)
            h2 = m * sqrtx / 2 / (l - x^2) *
                 ((l - x^2) * (n * pi / 2 + atan(sqrtx)) / sqrtx - (l + x))

            B = 27 * h2 / 4 / (sqrtx * (1 + h1))^3
            f = 0.0
            if B < 0.0
                f = 2 * cos(acos(sqrt(B + 1)) / 3)
            else
                A = (sqrt(B) + sqrt(B + 1))^(1.0 / 3.0)
                f = A + 1.0 / A
            end
            y = 2 / 3 * sqrtx * (1 + h1) * (sqrt(B + 1) / f + 1)

            xold = x
            x = 0.5 * ((m / y^2 - (1 + l)) - sqrt(max((m / y^2 - (1 + l))^2 - 4 * l, 0.0)))
            if x < 0
                return (zeros(3), zeros(3))
            end
            loops += 1
        end

        a = s * (1 + λ)^2 * (1 + x) * (l + x) / 8 / x

        p = 2 * r1n * r2n * sin(dθ / 2)^2 * ((1 + x) / (l + x)) / s / (1 + λ)^2

        e = sqrt(1 - p / a)

    else # Standard case
        x = n == 0 ? l : 1 + 4 * l
        xold = 10 * l
        y = 0.0
        loops = 0
        while abs(x - xold) > 1e-16 && loops < 50
            h1 = 0.0
            h2 = 0.0
            if n > 0
                temp = (n * pi / 2 + atan(sqrt(x))) / sqrt(x)
                temp2 = 4 * x^2 * (1 + 2 * x + l)

                h1 = (l + x)^2 / temp2 * (3 * (1 + x)^2 * temp - (3 + 5 * x))
                h2 = m / temp2 * ((x^2 - (1 + l) * x - 3 * l) * temp + (3 * l + x))
            else
                η = x / ((sqrt(1 + x) + 1)^2)
                a0 = [8 * (sqrt(1 + x) + 1), 1.0, 9 * η / 7]
                b0 = [0.0, 3.0, 5 + η]

                xi = continuedfraction(a0, b0, x -> a1(η, x), b1)

                h1 = (l + x)^2 * (1 + 3 * x + xi)
                h1 = h1 / (1 + 2 * x + l) / (4 * x + xi * (3 + x))

                h2 = m * (x - l + xi) / (1 + 2 * x + l) / (4 * x + xi * (3 + x))
            end

            # Alternate cubic solve
            B = 27.0 * h2 / 4 / (1 + h1)^3
            U = B / 2 / (sqrt(1 + B) + 1)

            a0 = [1 / 3, 4 / 27 * U]
            b0 = [0.0, 1.0]
            K = continuedfraction(a0, b0, x -> a2(U, x), b1)

            y = (1 + h1) / 3 * (2 + sqrt(1 + B) / (1 + 2 * U * K^2))

            xold = x
            x = sqrt(((1.0 - l) / 2)^2 + m / y / y) - (1.0 + l) / 2
            if x < 0.0
                return (zeros(3), zeros(3))
            end
            loops += 1
        end
        temp = sin(dθ / 2)^2
        p = 2 * r1n * r2n * y^2 * (1 + x)^2 * temp / m / s / (1 + λ)^2
        ϵ = (r2n - r1n) / r1n
        e = sqrt((ϵ^2 + 4 * r2n / r1n * temp * ((l - x) / (l + x))^2) /
                 (ϵ^2 + 4 * r2n / r1n * temp))
    end

    # Hodograph solution to resolve plane ambiguities
    A = mu * (1 / r1n - 1 / p)
    B = (mu * e / p)^2 - A^2

    x1 = B > 0 ? -sqrt(B) : 0.0

    node = zeros(3)
    if abs(sdθ) < 1e-8
        node = unit(cross(r1, v1))
        if e < 1
            Ptx = 2 * pi * sqrt(p^3 / mu / (1 - e^2)^3)
            if mod(tof, Ptx) > Ptx / 2
                x1 = -x1
            end
        end
    else
        node = unit(cross(r1, r2))
        if mod(dθ, 2 * pi) > pi
            node = -node
        end
        y2a = mu / p - x1 * sdθ + A * cdθ
        y2b = mu / p + x1 * sdθ + A * cdθ

        if abs(mu / r2n - y2b) < abs(mu / r2n - y2a)
            x1 = -x1
        end
    end

    v1new = sqrt(mu * p) / r1n * (x1 / mu * r1 + 1 / r1n * cross(node, r1))
    x2 = x1 * cdθ + A * sdθ
    v2 = sqrt(mu * p) / r2n * (x2 / mu * r2 + 1 / r2n * cross(node, r2))

    return (v1new, v2)
end
