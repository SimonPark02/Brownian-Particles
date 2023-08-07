"""
    This is a simulation for a single active brownian particle in 2D.
    The Langevin equation is numerically integrated by Euler method.
    Author: parkyongjun@snu.ac.kr
"""

using Distributions

const m = 1e-9
const R = 1e-7
const eta = 0.001
const kB = 1.38e-23
const T = 300
const gamma = 6 * pi * eta * R
const D = kB * T / gamma
const Dr = kB * T / (8 * pi * eta * R^3)
const V0 = 2e-6
const dt = 1e-3
const MAX_IT = 100000

function step!(X::Vector{Float64}, V::Vector{Float64}, A::Vector{Float64})
    for i=1:2
        A[i] = rand(Normal(0, 1))
    end
    X[1] += dt * V[1]
    X[2] += dt * V[2]
    X[3] += dt * sqrt(2 * Dr) * rand(Normal(0, 1))
    V[1] += dt * (- gamma * V[1] + gamma * V0 * cos(X[3]) + sqrt(2 * gamma * kB * T) * A[1]) / m
    V[2] += dt * (- gamma * V[2] + gamma * V0 * sin(X[3]) + sqrt(2 * gamma * kB * T) * A[2]) / m
end

function save2csv(f::IOStream, t::Float64, X::Vector{Float64}, V::Vector{Float64})
    write(f, string(t, ", "))
    for i=1:3
        write(f, string(X[i], ", "))
    end
    for i=1:2
        write(f, string(V[i], ", "))
    end
    write(f, "\n")
end

function main()
    X = [0.0, 0.0, 2 * pi * rand()]
    V = rand(Float64, 2) * 1e-6
    A = [0.0, 0.0]
    f = open(joinpath(pwd(), "single_active_bp.csv"), "w")
    for i=1:MAX_IT
        step!(X, V, A)
        save2csv(f, i * dt, X, V)
    end
    close(f)
end

@time main()