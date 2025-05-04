temp = load("src/results.jld2")
results = temp["results"]


##########################
# 1: 5 deg GEO transfer E
x0 = [-30726.027555, 28903.307, -3.75426, -2.10596, -2.23904, -0.001]
xf0 = [-33128.188504, 26115.380289, 2.26439871, -1.902702, -2.4141689, -0.00149266]


############################
# 2: 10 deg xfer W
x0 = [-30726.027555, 28903.307, -3.75426, -2.10596, -2.23904, -0.001]
xf0 = [-25240.2420737, 33799.701807, -17.4391459, -2.4630497, -1.838969, -0.00013836]



#############################
# 3: 25 deg xfer E
x0 = [-30726.027555, 28903.307, -3.75426, -2.10596, -2.23904, -0.001]
xf0 = [-40062.2849349, 13209.954533, 19.83209, -0.96162606, -2.919813, -0.00378514]


# Note: all three cases are coplanar.

tof = 60.0 * 60 * 52


r, v = universalkepler(SA[xf0[1:3]...], SA[xf0[4:6]...], tof, 3.986e5)
xf = [r; v]
x1 = copy(xf)
mu = 3.986e5
n = _get_default_half_revs(0.5 * (norm(x0[1:3]) + norm(x1[1:3])), tof, mu=mu)
verbose = true
vxfer, vstop = basic_lambert(x0[1:3], x0[4:6], x1[1:3], tof, n, mu=mu, verbose=true, v2=x1[4:6])
N = Int(tof / 3600) #default 20
tN, wN = _get_nodes_and_weights(N)

# Get the dv numbers from this block
nodes = primer_lambert(x0, x1, tof, verbose=true, driftSearch=false)

model2, states2, control2 = pseudospectral_continuous_lambert(x0, x1, tof, poly_order=N, verbose=true, constrain_u=true, mu=mu)
dvtotal2 = 0.0
for i in 1:N+1
    dvtotal2 += wN[i] * (norm(control2[i, :]))
end
dvtotal2 *= tof / 2
@info "DV total is $dvtotal2"

nodes2 = pseudospectral_impulsive_lambert(x0, x1, tof, poly_order=N, verbose=true, constrain_u=true)
dv = 0.0
for node in nodes2
    dv += norm(node.outVelocity - node.inVelocity)
end
@info "DV total is $dv"

