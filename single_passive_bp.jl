"""
    This is a simulation for a single passive brownian particle in 2D.
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
const dt = 1e-3
const MAX_IT = 100000

function step!(X::Vector{Float64}, V::Vector{Float64}, A::Vector{Float64})
    for i=1:2
        A[i] = rand(Normal(0, 1))
    end
    X[:] += dt * V[:]
    V[:] += dt * (- gamma * V[:] + sqrt(2 * gamma * kB * T) * A[:]) / m
end

function save2csv(f::IOStream, t::Float64, X::Vector{Float64}, V::Vector{Float64})
    write(f, string(t, ", "))
    for i=1:2
        write(f, string(X[i], ", "))
    end
    for i=1:2
        write(f, string(V[i], ", "))
    end
    write(f, "\n")
end

function main()
    X = [0.0, 0.0]
    V = rand(Float64, 2) * 1e-6
    A = [0.0, 0.0]
    f = open(joinpath(pwd(), "single_passive_bp.csv"), "w")
    for i=1:MAX_IT
        step!(X, V, A)
        save2csv(f, i * dt, X, V)
    end
    close(f)
end

@time main()