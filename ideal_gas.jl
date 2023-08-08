"""
    This is a simulation for 2D ideal gas confined in a rectangular geometry.
    Author: parkyongjun@snu.ac.kr
"""

const N = 2000
const V = 500.0
const dt = 8e-5
const R = 3e-2

function init()
    Gas = rand(Float64, (N, 4))
    for i in 1:N
        v = V * rand()
        theta = 2 * pi * rand()
        Gas[i, 3] = v * cos(theta)
        Gas[i, 4] = v * sin(theta)
    end
    return Gas
end

function go!(Gas::Matrix{Float64})
    for i = 1:2
        Gas[:, i] += Gas[:, i + 2] * dt
    end
end

function collision_wall!(Gas::Matrix{Float64})
    for j in 1:2
        for i in 1:N
            if (abs(Gas[i, j] - 0.5) > (1/2 - R)) && ((Gas[i, j] - 0.5) * Gas[i, j + 2] > 0)
                Gas[i, j + 2] = - Gas[i, j + 2]
            end
        end
    end
end

function collision_two!(Gas::Matrix{Float64}, d::Vector{Float64}, vr::Vector{Float64})
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

function step!(Gas::Matrix{Float64}, d::Vector{Float64}, vr::Vector{Float64})
    go!(Gas)
    collision_wall!(Gas)
    collision_two!(Gas, d, vr)
end

function save2csv(f::IOStream, t::Float64, Gas::Matrix{Float64})
    write(f, string(t, ", "))
    for i=1:N
        write(f, string(sqrt(Gas[i, 3] * Gas[i, 3] + Gas[i, 4] * Gas[i, 4]), ", "))
    end
    write(f, "\n")
end

function main()
    Gas = init()
    d = [0.0, 0.0]
    vr = [0.0, 0.0]
    f = open(joinpath(pwd(), "ideal_gas.csv"), "w")
    for i=1:1000
        step!(Gas, d, vr)
        save2csv(f, i * dt, Gas)
    end
    close(f)
end

@time main()
