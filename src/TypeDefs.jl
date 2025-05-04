export PrimerNode
mutable struct PrimerNode
    time::Float64
    position::SVector{3,Float64}
    inVelocity::SVector{3,Float64}
    outVelocity::SVector{3,Float64}
end

function Base.isless(a::PrimerNode, b::PrimerNode)
    return a.time < b.time
end

struct Particle
    position::Vector{Float64}
    velocity::Vector{Float64}
    bestPos::Vector{Float64}
    bestScore::Float64
end

export statePrimerHist
struct statePrimerHist
    t::Vector{Float64}
    state::Matrix{Float64}
    primers::Matrix{Float64}
end
