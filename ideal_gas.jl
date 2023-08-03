"""
    This is a simulation for 2D ideal gas confined in a rectangular geometry.
"""

function init(N::Int64, V::Float64)
    Gas = rand(Float64, (N, 4))
    for i in 1:N
        v = V * rand()
        theta = 2 * pi * rand()
        Gas[i, 3] = v * cos(theta)
        Gas[i, 4] = v * sin(theta)
    end
    return Gas
end

function go!(Gas::Matrix{Float64}, dt::Float64)
    for i = 1:2
        Gas[:, i] += Gas[:, i + 2] * dt
    end
end

function collision_wall!(Gas::Matrix{Float64}, N::Int64, R::Float64)
    for j in 1:2
        for i in 1:N
            if (abs(Gas[i, j] - 0.5) > (1/2 - R)) && ((Gas[i, j] - 0.5) * Gas[i, j + 2] > 0)
                Gas[i, j + 2] = - Gas[i, j + 2]
            end
        end
    end
end

function collision_two!(Gas::Matrix{Float64}, N::Int64, R::Float64, d::Vector{Float64}, vr::Vector{Float64})
    for i in 1:(N-1)
        for j in (i+1):N
            for k in 1:2
                d[k] = Gas[i, k] - Gas[j, k]
                vr[k] = Gas[i, k + 2] - Gas[j, k + 2]
            end
            dd = 0
            dvr = 0
            for k = 1:2
                dd += d[k] * d[k]
                dvr += d[k] * vr[k]
            end
            if (sqrt(dd) < 2 * R) && (dvr < 0)
                d = - dvr / dd * d
                for k = 1:2
                    Gas[i, k + 2] += d[k]
                    Gas[j, k + 2] -= d[k]
                end
            end
        end
    end
end

function step!(Gas::Matrix{Float64}, N::Int64, dt::Float64, R::Float64, d::Vector{Float64}, vr::Vector{Float64})
    go!(Gas, dt)
    collision_wall!(Gas, N, R)
    collision_two!(Gas, N, R, d, vr)
end

function save2csv(f::IOStream, t::Float64, Gas::Matrix{Float64}, N::Int64)
    write(f, string(t, ", "))
    for i=1:N
        write(f, string(sqrt(Gas[i, 3] * Gas[i, 3] + Gas[i, 4] * Gas[i, 4]), ", "))
    end
    write(f, "\n")
end

function main(N::Int64, V::Float64, dt::Float64, R::Float64)
    Gas = init(N, V)
    d = [0.0, 0.0]
    vr = [0.0, 0.0]
    f = open(joinpath(pwd(), "dist.csv"), "w")
    for i=1:1000
        step!(Gas, N, dt, R, d, vr)
        save2csv(f, i * dt, Gas, N)
    end
    close(f)
end

@time main(2000, 500.0, 8e-5, 3e-2)
