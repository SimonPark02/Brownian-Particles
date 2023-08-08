using Plots

const epsilon = 1.0
const sigma = 10.0

function U(r::Float64)
    return 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6)
end

function F(r::Float64)
    return 24 * epsilon * sigma / r^2 * (2 * (sigma / r)^11 - (sigma / r)^5) 
end

plot(F, 9, 15)