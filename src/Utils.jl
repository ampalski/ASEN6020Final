export unit
function unit(r1::AbstractVector)
    return (1 / norm(r1)) .* r1
end

function d8rt(x)
    return sqrt(sqrt(sqrt(x)))
end

export anglevec
function anglevec(r1::AbstractVector, r2::AbstractVector)
    return acos(dot(r1, r2) / norm(r1) / norm(r2))
end

"""
Solves a generalized continued fraction given starting a0 and b0 values and 
algorithms for later a and b values, using the Lentz algorithm

Must be in the form `value = b0 + a1 / (b1 + a2 / (b2 + ...`
"""
function continuedfraction(
    a0::AbstractVector,
    b0::AbstractVector,
    a::Function,
    b::Function,
)
    tiny = 1e-30

    #isolate b0
    b00 = popfirst!(b0)

    #Set up initial values
    f = b00 == 0 ? tiny : b00

    C = f
    D = 0

    ind = 1
    del = 100 # C * D, should approach 1 over time

    while abs(del - 1) > eps()
        bj = ind > length(b0) ? b(ind) : b0[ind]
        aj = ind > length(a0) ? a(ind) : a0[ind]

        D = bj + aj * D
        D = D == 0 ? tiny : D

        C = bj + aj / C
        C = C == 0 ? tiny : C

        D = 1 / D
        del = C * D
        f = del * f
        ind += 1
    end

    return f
end

export _get_nodes_and_weights
function _get_nodes_and_weights(N::Int)
    return gausslobatto(N + 1)
end

export _differentiation_matrix
function _differentiation_matrix(tN)
    N = length(tN) - 1
    kl0 = -N * (N + 1.0) / 4.0
    klN = N * (N + 1.0) / 4.0
    D = zeros(N + 1, N + 1)

    for k in 1:N+1
        for l in 1:N+1
            if k == l && k == 1
                D[k, l] = kl0
            elseif k == l && k == N + 1
                D[k, l] = klN
            elseif k != l
                D[k, l] = Pl(tN[k], N) / Pl(tN[l], N) / (tN[k] - tN[l])
            end
        end
    end
    return D
end
