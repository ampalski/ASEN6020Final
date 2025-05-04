function TLamb(m, q, qsqfm1, x, n)
    # Initialize variables
    sw = 0.4
    LM1 = n == -1
    L1 = n >= 1
    L2 = n >= 2
    L3 = n == 3
    qsq = q * q
    xsq = x * x
    U = (1 - x) * (1 + x)

    T = 0.0
    DT = 0.0
    D2T = 0.0
    D3T = 0.0

    # Direct computations (not series)
    if LM1 || m > 0 || x < 0 || abs(U) > sw
        y = sqrt(abs(U))
        z = sqrt(qsqfm1 + qsq * xsq)
        qx = q * x
        a = b = aa = bb = 0.0
        if qx <= 0
            a = z - qx
            b = q * z - x
        end
        if qx < 0 && LM1
            aa = qsqfm1 / a
            bb = qsqfm1 * (qsq * U - xsq) / b
        end
        if qx == 0 && LM1 || qx > 0
            aa = z + qx
            bb = q * z + x
        end
        if qx > 0
            a = qsqfm1 / aa
            b = qsqfm1 * (qsq * U - xsq) / bb
        end
        if !LM1
            if qx * U >= 0
                g = x * z + q * U
            else
                g = (xsq - qsq * U) / (x * z - q * U)
            end
            f = a * y
            if x <= 1
                T = m * pi + atan(f, g)
            else
                if f > sw
                    T = log(f + g)
                else
                    fg1 = f / (g + 1)
                    term = 2 * fg1
                    fg1sq = fg1 * fg1
                    T = term
                    ii = 1
                    Told = 2 * T
                    iter = 0
                    while T != Told && iter < 100 #this block should be atanh ... see if I can condense
                        ii += 2
                        term = term * fg1sq
                        Told = T
                        T += T + term / ii
                        iter += 1
                    end
                end
            end
            T = 2 * (T / y + b) / U
            if L1 && z != 0
                qz = q / z
                qz2 = qz * qz
                qz = qz * qz2

                DT = (3 * x * T - 4 * (a + qx * qsqfm1) / z) / U
                L2 && (D2T = (3 * T + 5 * x * DT + 4 * qz * qsqfm1) / U)
                L3 && (D3T = (8 * DT + 7 * x * D2T - 12 * qz * qz2 * x * qsqfm1) / U)
            end
        else
            DT = b
            D2T = bb
            D3T = aa
        end
    else # Compute by series
        u0 = 1.0
        L1 && (u1 = 1.0)
        L2 && (u2 = 1.0)
        L3 && (u3 = 1.0)
        term = 4.0
        tq = q * qsqfm1
        ii = 0
        tqsum = q < 0.5 ? 1 - q * qsq : (1 / (1 + q) + q) * qsqfm1
        ttmold = term / 3
        T = ttmold * tqsum
        Told = T * 2
        while ii < n || T != Told
            ii += 1
            p = ii
            u0 = u0 * ii
            L1 && ii > 1 && (u1 = u1 * U)
            L2 && ii > 2 && (u2 = u2 * U)
            L3 && ii > 3 && (u3 = u3 * U)
            term = term * (p - 0.5) / p
            tq = tq * qsq
            tqsum += tq
            Told = T
            tterm = term / (2 * p + 3)
            tqterm = tterm * tqsum
            T = T - u0 * ((1.5 * p + 0.25) * tqterm / (p * p - 0.25) - ttmold * tq)
            ttmold = tterm
            tqterm = tqterm * p

            L1 && (DT = DT + tqterm * u1)
            L2 && (D2T = D2T + tqterm * u2 * (p - 1))
            L3 && (D3T = D3T + tqterm * u3 * (p - 1) * (p - 2))
        end
        #@show ii
        L3 && (D3T = 8 * x * (1.5 * D2T - xsq * D3T))
        L2 && (D2T = 2 * (2 * xsq * D2T - DT))
        L1 && (DT = -2 * x * DT)
        T = T / xsq
    end
    return (T, DT, D2T, D3T)
end

function XLamb(m, q, qsqfm1, Tin)
    #Initialize variables
    tol = 3e-7
    C0 = 1.7
    C1 = 0.5
    C2 = 0.03
    C3 = 0.15
    C41 = 1.0
    C42 = 0.24
    thr2 = atan(qsqfm1, 2 * q) / pi
    check2 = false

    #single revolution starter from T (at x=0) & Bilinear
    if m == 0
        n = 1
        T0, DT, D2T, D3T = TLamb(m, q, qsqfm1, 0.0, 0.0)
        TDiff = Tin - T0
        if TDiff <= 0
            x = T0 * TDiff / (-4 * Tin) #-4 is the value for DT for x=0
        else
            x = -TDiff / (TDiff + 4)
            w = x + C0 * sqrt(2 * (1 - thr2))
            if w < 0
                x = x - sqrt(d8rt(-w)) * (x + sqrt(TDiff / (TDiff + 1.5 * T0)))
            end
            w = w = 4 / (4 + TDiff)
            x = x * (1.0 + x * (C1 * w - C2 * x * sqrt(w)))
        end
        xpl = x
    else #Multi-rev starters
        xm = 1 / (1.5 * (m + 0.5) * pi)
        xm = thr2 < 0.5 ? d8rt(2 * thr2) * xm : (2 - d8rt(2 - 2 * thr2)) * xm

        check = false
        Tmin = 0
        for ii in 1:12
            Tmin, DT, D2T, D3T = TLamb(m, q, qsqfm1, xm, 3)
            if D2T == 0
                check = true
                xtest = 0
            else
                xmold = xm
                xm -= DT * D2T / (D2T * D2T - DT * D3T / 2)
                xtest = abs(xmold / xm - 1.0)
            end
            if xtest < tol || check
                check = true
                break
            end
        end
        #@info ii

        if !check #if Tmin is never found
            n = -1
            return (n, 0.0, 0.0)
        end

        TDiffm = Tin - Tmin
        if TDiffm < 0 #no solution with this m value
            n = 0
            return (n, 0.0, 0.0)
        elseif TDiffm == 0 #unique solution
            x = xm
            xpl = x
            n = 1
            return (n, x, xpl)
        else
            n = 3
            D2T == 0 && (D2T = 6 * m * pi)
            x = sqrt(TDiffm / (D2T / 2 + TDiffm / (1 - xm)^2))
            w = xm + x
            w = w * 4 / (4 + TDiffm) + (1 - w)^2
            x = x * (1 - (1 + m + C41 * (thr2 - 0.5)) / (1 + C3 * m) *
                         x * (C1 * w + C2 * x * sqrt(w))) + xm
            D2T2 = D2T / 2
            if x >= 1
                n = 1
                check2 = true
            end
        end

    end
    # Have a starter, proceed by Halley
    jj = 0
    while jj < 10 #prevent stuck in "GO TO 5" loop
        if !check2
            for _ in 1:3
                T, DT, D2T, D3T = TLamb(m, q, qsqfm1, x, 2)
                T = Tin - T
                DT != 0.0 && (x += T * DT / (DT * DT + T * D2T / 2))
            end
            if n != 3
                return (n, x, xpl)
            end
            n = 2
            xpl = x
        end

        T0, DT, D2T, D3T = TLamb(m, q, qsqfm1, 0.0, 0)
        TDiff0 = T0 - Tmin
        TDiff = Tin - T0
        if TDiff <= 0
            x = xm - sqrt(TDiffm / (D2T2 - TDiffm * (D2T2 / TDiff0 - 1 / xm^2)))
        else
            x = -TDiff / (TDiff + 4)
            w = x + C0 * sqrt(2 * (1 - thr2))
            if w < 0
                x = x - sqrt(d8rt(-w)) * (x + sqrt(TDiff / (TDiff + 1.5 * T0)))
            end
            w = 4 / (4 + TDiff)
            x = x * (1 + (1 + m + C42 * (thr2 - 0.5)) / (1 + C3 * m) *
                         x * (C1 * w - C2 * x * sqrt(w)))
            if x <= -1
                n = n - 1
                n == 1 && (x = xpl)
            end
        end
        jj = jj + 1
        check2 = false
    end
    @show jj
    error("Too many iterations in VLamb")
end

function VLamb(mu, r1, r2, theta, Tdelt)
    #Initialize Vars
    VR = zeros(2, 2)
    VT = zeros(2, 2)
    twopi = 2 * pi
    m = floor(theta / twopi)
    thr2 = theta / 2 - m * pi
    dr = r1 - r2
    r1r2 = r1 * r2
    r1r2th = 4 * r1r2 * sin(thr2)^2
    csq = dr^2 + r1r2th
    c = sqrt(csq)
    s = (r1 + r2 + c) / 2

    gms = sqrt(mu * s / 2)
    qsqfm1 = c / s
    q = sqrt(r1r2) * cos(thr2) / s

    if c != 0
        rho = dr / c
        sig = r1r2th / csq
    else
        rho = 0.0
        sig = 1.0
    end

    # Find x
    T = 4 * gms * Tdelt / s / s
    n, x1, x2 = XLamb(m, q, qsqfm1, T)

    # Solve for velocities
    for ii in 1:n
        x = ii == 1 ? x1 : x2

        _, qzminx, qzplx, zplqx = TLamb(m, q, qsqfm1, x, -1)
        vt2 = gms * zplqx * sqrt(sig)
        vr1 = gms * (qzminx - qzplx * rho) / r1
        vt1 = vt2 / r1
        vr2 = -gms * (qzminx + qzplx * rho) / r2
        vt2 = vt2 / r2

        VR[ii, 1] = vr1
        VT[ii, 1] = vt1
        VR[ii, 2] = vr2
        VT[ii, 2] = vt2
    end
    return (VR, VT)
end

# Manual translation from Gooding's (1989) original Fortran
# Replacements of various GOTO statements may need work
function lambert_gooding(
    r1::AbstractVector, #Starting position vector
    r2::AbstractVector, #Ending position vector
    tof::AbstractFloat, #time of flight
    n::Int; # Number of half revolutions between r1 and r2
    mu::AbstractFloat=3.986e5, #grav parameter for units used
    pick::Int=0 #which velocity to return (1 or 2), or 0 for both
)
    theta = anglevec(r1, r2)
    if iseven(n)
        theta = theta + n * pi
    else
        theta = (n + 1) * pi - theta
    end

    VR, VT = VLamb(mu, norm(r1), norm(r2), theta, tof)

    # Rotate velocities from radial/transverse to ECI

    if iseven(n)
        C = unit(cross(r1, r2))
    else
        C = unit(cross(r2, r1))
    end

    R = unit(r1)
    I = cross(C, R)
    v11 = [R I C] * [VR[1, 1], VT[1, 1], 0.0]
    v12 = [R I C] * [VR[2, 1], VT[2, 1], 0.0]

    R = unit(r2)
    I = cross(C, R)
    v21 = [R I C] * [VR[1, 2], VT[1, 2], 0.0]
    v22 = [R I C] * [VR[2, 2], VT[2, 2], 0.0]

    if pick == 1
        return (v11, v21)
    elseif pick == 2
        return (v12, v22)
    else
        return ([v11, v12], [v21, v22])
    end
end
