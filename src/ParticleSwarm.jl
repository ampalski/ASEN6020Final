function InitSwarm(bounds, numParticles::Int, evalFun)
    swarm = Vector{Particle}(undef, numParticles)

    for i in 1:numParticles
        position = Vector{Float64}(undef, length(bounds))
        velocity = copy(position)

        #Uniform initialization of position and velocity
        for j in eachindex(bounds)
            r = bounds[j][2] - bounds[j][1]
            position[j] = bounds[j][1] + r * rand()
            velocity[j] = -r + 2.0 * r * rand()
        end

        #Initial Evaluation
        score = evalFun(position)

        #Save off
        swarm[i] = Particle(position, velocity, position, score)
    end

    return swarm
end

function findBest(swarm::Vector{Particle})
    bestPosition = swarm[1].bestPos
    bestScore = swarm[1].bestScore

    for i in 2:length(swarm)
        score = swarm[i].bestScore
        if score < bestScore
            bestScore = score
            bestPosition = swarm[i].bestPos
        end
    end
    return (bestPosition, bestScore)
end

function IterateSwarm(
    swarm::Vector{Particle},
    bounds::Vector{Vector{Float64}},
    W::Float64,
    y1::Float64,
    y2::Float64,
    numNeighbors::Int,
    evalFun,
    gradientFrac::Float64,
    gradientFun,
)

    newSwarm = Vector{Particle}(undef, length(swarm))
    for i = 1:length(swarm)
        neighbors = Vector{Particle}(undef, numNeighbors)
        inds = []
        particle = swarm[i]

        while length(inds) < numNeighbors
            randInd = rand(1:length(swarm))
            if randInd == i || in(randInd, inds)
                continue
            end

            push!(inds, randInd)

            neighbors[length(inds)] = swarm[randInd]
        end


        # Find the best score of the neighbors
        (bestNeighbor, bestNeighborScore) = findBest(neighbors)

        # Update the velocity
        v = zeros(length(bounds))

        # Experimental
        if rand() < gradientFrac
            v = -rand(length(bounds)) .* gradientFun(particle.position)
        else
            u1 = rand(length(bounds))
            u2 = rand(length(bounds))
            v = W .* particle.velocity + y1 .* u1 .* (particle.bestPos - particle.position) +
                y2 .* u2 .* (bestNeighbor - particle.position)
        end

        # Update the position
        x = particle.position + v

        # Enforce bounds
        for j in 1:length(bounds)
            if x[j] > bounds[j][2]
                x[j] = bounds[j][2]
                if v[j] > 0
                    v[j] = 0
                end
            end

            if x[j] < bounds[j][1]
                x[j] = bounds[j][1]
                if v[j] < 0
                    v[j] = 0
                end
            end
        end

        # Evaluate Particle
        score = evalFun(x)

        # Check if this is the best this particle has seen
        if score < particle.bestScore
            newSwarm[i] = Particle(x, v, x, score)
        else
            newSwarm[i] = Particle(x, v, particle.bestPos, particle.bestScore)
        end


    end
    return newSwarm
end

export ParticleSwarmOptimizer
"""
    bestAnswer = ParticleSwarmOptimizer(bounds, evalFun; <keyword arguments>)

Optimize (minimize) a function using a Particle Swarm.

# Arguments
- `bounds::AbstractVector{AbstractVector{Float64}}`: An n-length vector of 2-length vectors
that describe the lower and upper bounds for each element of the solution space
- `evalFun`: A function that accepts an n-length vector as an argument that evaluates the
function to be minimized
- `numParticles::Int = 100`: The number of particles considered during each iteration
- `maxIter::Int = 200`: The maximum iterations taken during optimation. Note that there are
not currently any short-circuit methods to end optimization early
- `numNeighborsFrac::Float64 = 0.25`: The fraction of the full population numParticles
considered by default to be part of a given particle's neighborhood
- `selfWeight::Float64 = 1.49`: The weight given to the velocity update based on the best
value found by that individual particle
- `socialWeight::Float64 = 1.49`: The weight given to the velocity update based on the best
value found by any particle in the neighborhood
- `inertiaRange::AbstractVector{Float64} = [0.1, 1.1]`: A 2-vector bounding the weight given
to the velocity update based on the particle's previous velocity
- `verbose::Bool = false`: A flag to tell the optimizer to output status to the REPL
"""
function ParticleSwarmOptimizer(
    bounds::Vector{Vector{Float64}},
    evalFun;
    numParticles::Int=100,
    maxIter::Int=200,
    numNeighborsFrac::Float64=0.25,
    selfWeight::Float64=1.49,
    socialWeight::Float64=1.49,
    inertiaRange::AbstractVector{Float64}=[0.1, 1.1],
    verbose::Bool=false,
    logging::Bool=false,
    gradientFrac::Float64=0.0,
    gradientFun=sin,
)
    if logging
        logfile = FormatLogger(open("log.txt", "w")) do io, args
            # Write the module, level and message only
            println(io, args._module, " | ", "[", args.level, "] ", args.message)
        end
    end
    #Initialize swarm
    swarm = InitSwarm(bounds, numParticles, evalFun)

    #Initialize helper vars and weights
    N = max(2, floor(Int, numParticles * numNeighborsFrac)) #neighborhood size
    minNeighborhoodSize = N
    W = maximum(inertiaRange)
    c = 0
    y1 = selfWeight
    y2 = socialWeight
    (bestPosition, bestScore) = findBest(swarm)

    #Iterate
    for i in 1:maxIter
        if verbose
            println("Running Iteration $i of $maxIter")
        end
        swarm = IterateSwarm(swarm, bounds, W, y1, y2, N, evalFun, gradientFrac, gradientFun)

        # Check if this is the best the swarm has seen
        flag = false
        (checkBestPosition, checkBestScore) = findBest(swarm)
        #println("Best score $checkBestScore at $checkBestPosition")
        if checkBestScore < bestScore
            bestPosition = checkBestPosition
            bestScore = checkBestScore
            flag = true

            verbose && display(bestScore)
        end

        # Update neighborhood parameters
        if flag
            c = max(0, c - 1)
            N = minNeighborhoodSize
            if c < 2
                W *= 2.0
            elseif c > 5
                W /= 2.0
            end

            W = clamp(W, inertiaRange[1], inertiaRange[2])

        else
            c += 1
            N = min(N + minNeighborhoodSize, numParticles - 1)
        end

        if logging
            j = 0
            for p in swarm
                j += 1
                with_logger(logfile) do
                    @debug "Iter $i Particle $j Best Score: $(p.bestScore)"
                end
            end
        end
    end

    return bestPosition
end
