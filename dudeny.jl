module Dudeny

using Plots, JuMP, Ipopt

function xy_to_polar(x, y)
    θ = atan(y, x)
    r = hypot(x, y)
    return θ, r
end

function polar_to_xy(θ, r)
    return r * cos(θ), r * sin(θ)
end

function deg_to_rad(d)
    return 2π * (d / 360.0)
end

function equitriangle(side)
    xs = Vector{Float64}(undef, 4)
    ys = Vector{Float64}(undef, 4)
    xs[1] = 0.
    ys[1] = 0.
    xs[2] = xs[1] + side
    ys[2] = 0.
    x, y = polar_to_xy(deg_to_rad(120), side) .+ (side, 0)
    xs[3] = x
    ys[3] = y
    xs[4] = 0.
    ys[4] = 0.
    return xs, ys
end

function mid(xs, ys)
    return sum(xs) / 2, sum(ys) / 2
end

function circle(radius, c=(0., 0.); n_points=1000)
    θ = range(start=0, stop=2π, length=n_points)
    r = fill(radius, n_points)
    xy = polar_to_xy.(θ, r)
    xs = [p[1] + c[1] for p in xy]
    ys = [p[2] + c[2] for p in xy]
    return xs, ys
end

function arc(radius, c=(0., 0.); p1, p2, n_points=500, adj=0.0)
    p1 = p1 .- c
    p2 = p2 .- c
    θ₁, _ = xy_to_polar(p1...)
    θ₂, _ = xy_to_polar(p2...)
    θ = range(start=θ₁, stop=θ₂, length=n_points) .+ adj
    r = fill(radius, n_points)
    xy = polar_to_xy.(θ, r)
    xs = [p[1] + c[1] for p in xy]
    ys = [p[2] + c[2] for p in xy]
    return xs, ys
end

function distance(xs , ys)
    return hypot(xs[1]-xs[2], ys[1]-ys[2])
end

function gradient(xs, ys)
    m = (ys[1] - ys[2]) / (xs[1] - xs[2])
    c = ys[1] - m * xs[1]
    return m, c
end

function y_intercept(m, (px, py))
    return py - m * px
end

function compute_f(m, c, ex, ey, r)
    model = Model(Ipopt.Optimizer)
    @variable(model, x >= ex)
    @variable(model, y >= ey)
    @NLconstraint(model, (x - ex)^2 + (y - ey)^2 == r ^ 2)
    @constraint(model, y == m * x + c)
    optimize!(model)
    return value(x), value(y)
end

function compute_h(m, c, gx, gy, r)
    model = Model(Ipopt.Optimizer)
    @variable(model, x >= 0)
    @variable(model, y >= gy)
    @NLconstraint(model, (x - gx)^2 + (y - gy)^2 == r ^ 2)
    @constraint(model, y == m * x + c)
    optimize!(model)
    return value(x), value(y)
end

function compute_j(m, c, ex, ey, r)
    model = Model(Ipopt.Optimizer)
    @variable(model, x >= 0)
    @variable(model, y >= 0)
    @NLconstraint(model, (x - ex)^2 + (y - ey)^2 == r ^ 2)
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

function main(side::Number)
    xs, ys = equitriangle(side)
    dx, dy = mid(xs[3:4], ys[3:4])
    ex, ey = mid(xs[2:3], ys[2:3])

    m_ae, c_ae = gradient([ex,xs[1]], [ey,ys[1]])

    fx, fy = compute_f(m_ae, c_ae, ex, ey, side / 2.0)
    xs_ae = collect(range(start=0, stop=fx, length=500))
    ys_ae = m_ae .* xs_ae

    xs_bf, ys_bf = arc(side/ 2.0, (ex, ey), p1=(xs[3], ys[3]), p2=(fx, fy))

    gx, gy = mid([xs[1], fx], [ys[1], fy])

    m_eb, c_eb = gradient([ex, xs[3]], [ey, ys[3]])
    len_ag = distance([xs[1], gx], [ys[1], gy])
    xs_af, ys_af = arc(len_ag, (gx, gy), p1=(fx, fy), p2=(xs[1], ys[1]), adj=π)

    hx, hy = compute_h(m_eb, c_eb, gx, gy, len_ag)
    len_eh = distance([ex, hx], [ey, hy])
    jx, jy = compute_j(0, 0, ex, ey, len_eh)
    # xs_hj, ys_hj = arc(len_eh, (ex, ey), p1=(hx, hy), p2=(jx, jy), adj=π)
    kx, ky = jx + side / 2.0, 0.0

    m_je, c_je = gradient([jx, ex], [jy, ey])
    mp_je = -1.0 / m_je
    cp_je_d = y_intercept(mp_je, (dx, dy))
    lx, ly = lines_intersect(m_je, c_je, mp_je, cp_je_d)
    cp_je_k = y_intercept(mp_je, (kx, ky))
    mx, my = lines_intersect(m_je, c_je, mp_je, cp_je_k)

    plot(xs, ys, legend=false, aspect_ratio=:equal)
    plot!([lx, ex, xs[3], dx, lx],[ly, ey, ys[3], dy, ly], fill=(0, :lightgreen))
    plot!([xs[1], jx, lx, dx, xs[1]], [ys[1], jy, ly, dy, ys[1]], fill=(0, :red))
    plot!([jx, kx, mx, jx], [jy, ky, my, jy], fill=(0, :blue))
    plot!([kx, xs[2], ex, mx, kx], [ky, ys[2], ey, my, ky], fill=(0, :yellow))
    plot!([dx], [dy], marker=:circle)
    plot!([ex], [ey], marker=:circle)
    plot!([kx], [ky], marker=:circle)
    # plot!([fx], [fy], marker=:circle)
    # plot!([gx], [gy], marker=:circle)
    # plot!([hx], [hy], marker=:circle)
    # plot!([jx], [jy], marker=:circle)
     # plot!(xs_ae, ys_ae)
    # plot!(xs_bf, ys_bf)
    # plot!(xs_af, ys_af)
    # plot!(xs_hj, ys_hj)
    # plot!([xs[3], hx], [ys[3], hy])
    # plot!([jx, ex], [jy, ey])
    # plot!([dx, lx], [dy, ly])
    # plot!([kx, mx], [ky, my])

end

end # module
