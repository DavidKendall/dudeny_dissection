module Dudeny

using Plots, JuMP, Ipopt

function xy_to_polar((x, y))
    θ = atan(y, x)
    r = hypot(x, y)
    return θ, r
end

function polar_to_xy((θ, r))
    return r * cos(θ), r * sin(θ)
end

function deg_to_rad(d)
    return 2π * (d / 360.0)
end

function mid((ax,ay), (bx,by))
    return (ax + bx) / 2, (ay + by) / 2
end

function distance((ax, ay), (bx, by))
    return hypot(ax - bx, ay - by)
end

function gradient((ax, ay), (bx, by))
    return (ay - by) / (ax - bx)
end

function y_intercept(m, (px, py))
    return py - m * px
end

function equitriangle(side)
    a = 0.0, 0.0
    b = polar_to_xy((deg_to_rad(60), side))
    c = side, 0.0
    return (a, b, c)
end

function circle(radius, c=(0., 0.); n_points=1000)
    θ = range(start=0, stop=2π, length=n_points)
    r = fill(radius, n_points)
    points = zip(θ, r)
    xy = polar_to_xy.(points)
    xs = [p[1] + c[1] for p in xy]
    ys = [p[2] + c[2] for p in xy]
    return xs, ys
end

function arc(radius, c=(0., 0.); p1, p2, n_points=500, adj=0.0)
    p1 = p1 .- c
    p2 = p2 .- c
    θ₁, _ = xy_to_polar(p1)
    θ₂, _ = xy_to_polar(p2)
    θ = range(start=θ₁, stop=θ₂, length=n_points) .+ adj
    r = fill(radius, n_points)
    points = zip(θ, r)
    xy = polar_to_xy.(points)
    xs = [p[1] + c[1] for p in xy]
    ys = [p[2] + c[2] for p in xy]
    return xs, ys
end

function arc_line_intersect(m, c, (px, py), r; x_lower=0, y_lower=0)
    model = Model(Ipopt.Optimizer)
    @variable(model, x >= x_lower)
    @variable(model, y >= y_lower)
    @NLconstraint(model, (x - px)^2 + (y - py)^2 == r ^ 2)
    @constraint(model, y == m * x + c)
    optimize!(model)
    return value(x), value(y)
end

function lines_intersect(m1, c1, m2, c2)
    model = Model(Ipopt.Optimizer)
    @variable(model, x)
    @variable(model, y)
    @constraint(model, y == m1 * x + c1)
    @constraint(model, y == m2 * x + c2)
    optimize!(model)
    return value(x), value(y)
end

function rotate(p::Tuple{Real, Real}, c, α)
    p = p .- c
    θₚ, rₚ = xy_to_polar(p)
    θₚ += α
    return polar_to_xy((θₚ, rₚ)) .+ c
end

function rotate(v, c, α)
    return [rotate(p, c, α) for p in v]
end

function close(points)
    n = length(points) + 1
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    for i in 1:(n - 1)
        xs[i], ys[i] = points[i]
    end
    xs[n], ys[n] = points[1]
    return xs, ys
end

function compute_dissection(side)
    triangle = (a, b, c) = ((ax, ay), (bx, by), (cx, cy)) = equitriangle(side)
    d = mid(a, b)
    e = ex, ey = mid(b, c)
    m_ae = gradient(a, e)
    c_ae = y_intercept(m_ae, e)
    f = arc_line_intersect(m_ae, c_ae, e, side / 2.0, x_lower=ex, y_lower=ey)
    g = _, gy = mid(a, f)
    m_eb = gradient(e, b)
    c_eb = y_intercept(m_eb, e)
    len_ag = distance(a, g)
    h = arc_line_intersect(m_eb, c_eb, g, len_ag, y_lower=gy)
    len_eh = distance(e, h)
    j = jx, _ = arc_line_intersect(0, 0, e, len_eh)
    k = jx + side / 2.0, 0.0
    m_je = gradient(j, e)
    c_je = y_intercept(m_je, j)
    mp_je = -1.0 / m_je
    cp_je_d = y_intercept(mp_je, d)
    l = lines_intersect(m_je, c_je, mp_je, cp_je_d)
    cp_je_k = y_intercept(mp_je, k)
    m = lines_intersect(m_je, c_je, mp_je, cp_je_k)
    red = [a, j, l, d]
    green = [l, e, b, d]
    blue = [j, k, m]
    yellow = [k, c, e, m]
    return triangle, red, green, blue, yellow
end

function show_dissection(triangle, red, green, blue, yellow)
    xs, ys = close(triangle) 
    r_xs, r_ys = close(red) 
    b_xs, b_ys = close(blue) 
    g_xs, g_ys = close(green) 
    y_xs, y_ys = close(yellow) 
    (a, b, c) = triangle
    dx, dy = mid(a, b)
    ex, ey = mid(b, c)
    kx, ky = yellow[1]
    side = distance(a, b)
    limits = (-1, side+1)
    tick = floor(limits[1]):ceil(limits[2])
    plot(xs, ys, legend=false, aspect_ratio=:equal, tick=tick, limits=limits)
    plot!(r_xs, r_ys, fill=(0, :red))
    plot!(b_xs, b_ys, fill=(0, :blue))
    plot!(g_xs, g_ys, fill=(0, :lightgreen))
    plot!(y_xs, y_ys, fill=(0, :yellow))
    plot!([dx], [dy], marker=:circle)
    plot!([ex], [ey], marker=:circle)
    plot!([kx], [ky], marker=:circle)
end

function animate_dissection(triangle, red, green, blue, yellow)
    (a, b, c) = triangle
    d = dx, dy = mid(a, b)
    e = ex, ey = mid(b, c)
    k = kx, ky = yellow[1]
    xs, ys = close(triangle) 
    g_xs, g_ys = close(green) 
    α = 0.01
    total_angle = 0.0
    side = distance(a, b)
    limits = (-side/1.8, 2*side)
    tick = floor(limits[1]):ceil(limits[2])
    @gif while total_angle ≤ π    
        plot(xs, ys, legend=false, aspect_ratio=:equal, tick=tick, limits=limits)
        plot!(g_xs, g_ys, fill=(0, :lightgreen))
        red = rotate(red, d, -α)
        r_xs, r_ys = close(red) 
        plot!(r_xs, r_ys, fill=(0, :red))
        blue = rotate(blue, k, 2α)
        yellow = rotate(yellow, e, α)
        y_xs, y_ys = close(yellow) 
        plot!(y_xs, y_ys, fill=(0, :yellow))
        k = kx, ky = yellow[1]
        shift = k .- blue[2]
        blue = [p .+ shift for p in blue]
        b_xs, b_ys = close(blue) 
        plot!(b_xs, b_ys, fill=(0, :blue))
        plot!([dx], [dy], marker=:circle)
        plot!([ex], [ey], marker=:circle)
        plot!([kx], [ky], marker=:circle)
        total_angle += α
    end
end

function main(side)
    triangle, red, green, blue, yellow = compute_dissection(side)
    animate_dissection(triangle, red, green, blue, yellow)
    show_dissection(triangle, red, green, blue, yellow)
end

end # module
