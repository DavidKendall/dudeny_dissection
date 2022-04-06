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

function mid(xs, ys)
    return sum(xs) / 2, sum(ys) / 2
end

function distance(xs, ys)
    return hypot(xs[1]-xs[2], ys[1]-ys[2])
end

function gradient(xs, ys)
    return (ys[1] - ys[2]) / (xs[1] - xs[2])
end

function y_intercept(m, (px, py))
    return py - m * px
end

function equitriangle(side)
    xs = Vector{Float64}(undef, 4)
    ys = Vector{Float64}(undef, 4)
    xs[1] = 0.
    ys[1] = 0.
    xs[2] = xs[1] + side
    ys[2] = 0.
    x, y = polar_to_xy((deg_to_rad(120), side)) .+ (side, 0)
    xs[3] = x
    ys[3] = y
    xs[4] = 0.
    ys[4] = 0.
    return xs, ys
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

function compute_dissection(side)
    xs, ys = equitriangle(side)
    triangle = zip(xs, ys)
    dx, dy = mid(xs[3:4], ys[3:4])
    ex, ey = mid(xs[2:3], ys[2:3])
    m_ae = gradient([ex,xs[1]], [ey,ys[1]])
    c_ae = y_intercept(m_ae, (ex, ey))
    fx, fy = arc_line_intersect(m_ae, c_ae, (ex,ey), side / 2.0, x_lower=ex, y_lower=ey)
    gx, gy = mid([xs[1], fx], [ys[1], fy])
    m_eb = gradient([ex, xs[3]], [ey, ys[3]])
    c_eb = y_intercept(m_eb, (ex, ey))
    len_ag = distance([xs[1], gx], [ys[1], gy])
    hx, hy = arc_line_intersect(m_eb, c_eb, (gx, gy), len_ag, y_lower=gy)
    len_eh = distance([ex, hx], [ey, hy])
    jx, jy = arc_line_intersect(0, 0, (ex, ey), len_eh)
    kx, ky = jx + side / 2.0, 0.0
    m_je = gradient([jx, ex], [jy, ey])
    c_je = y_intercept(m_je, (jx, jy))
    mp_je = -1.0 / m_je
    cp_je_d = y_intercept(mp_je, (dx, dy))
    lx, ly = lines_intersect(m_je, c_je, mp_je, cp_je_d)
    cp_je_k = y_intercept(mp_je, (kx, ky))
    mx, my = lines_intersect(m_je, c_je, mp_je, cp_je_k)
    green = [(lx,ly), (ex,ey), (xs[3],ys[3]), (dx,dy), (lx, ly)]
    red = [(xs[1],ys[1]), (jx,jy), (lx,ly), (dx, dy), (xs[1], ys[1])]
    yellow = [(kx,ky), (xs[2],ys[2]), (ex,ey), (mx, my), (kx,ky)]
    blue = [(jx,jy), (kx,ky), (mx,my), (jx, jy)]
    return triangle, red, green, blue, yellow
end

function show_dissection(triangle, red, green, blue, yellow)
    xs, ys = [p[1] for p in triangle], [p[2] for p in triangle] 
    r_xs, r_ys = [p[1] for p in red], [p[2] for p in red]
    b_xs, b_ys = [p[1] for p in blue], [p[2] for p in blue]
    g_xs, g_ys = [p[1] for p in green], [p[2] for p in green]
    y_xs, y_ys = [p[1] for p in yellow], [p[2] for p in yellow]
    dx, dy = mid(xs[3:4], ys[3:4])
    ex, ey = mid(xs[2:3], ys[2:3])
    kx, ky = yellow[1]
    side = distance(xs[1:2], ys[1:2])
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
    xs, ys = [p[1] for p in triangle], [p[2] for p in triangle]
    dx, dy = mid(xs[3:4], ys[3:4])
    ex, ey = mid(xs[2:3], ys[2:3])
    g_xs, g_ys = [p[1] for p in green], [p[2] for p in green]
    α = 0.01
    total_angle = 0.0
    side = distance(xs[1:2], ys[1:2])
    limits = (-side/1.8, 2*side)
    tick = floor(limits[1]):ceil(limits[2])
    @gif while total_angle ≤ π    
        plot(xs, ys, legend=false, aspect_ratio=:equal, tick=tick, limits=limits)
        plot!(g_xs, g_ys, fill=(0, :lightgreen))
        red = rotate(red, (dx, dy), -α)
        r_xs, r_ys = [p[1] for p in red], [p[2] for p in red]
        plot!(r_xs, r_ys, fill=(0, :red))
        blue = rotate(blue, yellow[1], 2α)
        yellow = rotate(yellow, (ex, ey), α)
        y_xs, y_ys = [p[1] for p in yellow], [p[2] for p in yellow]
        plot!(y_xs, y_ys, fill=(0, :yellow))
        shift = (y_xs[1], y_ys[1]) .- blue[2]
        blue = [p .+ shift for p in blue]
        b_xs, b_ys = [p[1] for p in blue], [p[2] for p in blue]
        plot!(b_xs, b_ys, fill=(0, :blue))
        plot!([dx], [dy], marker=:circle)
        plot!([ex], [ey], marker=:circle)
        plot!([y_xs[1]], [y_ys[1]], marker=:circle)
        total_angle += α
    end
end

function main(side)
    triangle, red, green, blue, yellow = compute_dissection(side)
    animate_dissection(triangle, red, green, blue, yellow)
    show_dissection(triangle, red, green, blue, yellow)
end

# function main(side::Number)
#     xs, ys = equitriangle(side)
#     dx, dy = mid(xs[3:4], ys[3:4])
#     ex, ey = mid(xs[2:3], ys[2:3])

#     m_ae = gradient([ex,xs[1]], [ey,ys[1]])
#     c_ae = y_intercept(m_ae, (ex, ey))
#     fx, fy = arc_line_intersect(m_ae, c_ae, (ex,ey), side / 2.0, x_lower=ex, y_lower=ey)
#     xs_ae = collect(range(start=0, stop=fx, length=500))
#     ys_ae = m_ae .* xs_ae

#     xs_bf, ys_bf = arc(side/ 2.0, (ex, ey), p1=(xs[3], ys[3]), p2=(fx, fy))

#     gx, gy = mid([xs[1], fx], [ys[1], fy])

#     m_eb = gradient([ex, xs[3]], [ey, ys[3]])
#     c_eb = y_intercept(m_eb, (ex, ey))
#     len_ag = distance([xs[1], gx], [ys[1], gy])
#     xs_af, ys_af = arc(len_ag, (gx, gy), p1=(fx, fy), p2=(xs[1], ys[1]), adj=π)

#     hx, hy = arc_line_intersect(m_eb, c_eb, (gx, gy), len_ag, y_lower=gy)
#     len_eh = distance([ex, hx], [ey, hy])
#     jx, jy = arc_line_intersect(0, 0, (ex, ey), len_eh)
#     # xs_hj, ys_hj = arc(len_eh, (ex, ey), p1=(hx, hy), p2=(jx, jy), adj=π)
#     kx, ky = jx + side / 2.0, 0.0

#     m_je = gradient([jx, ex], [jy, ey])
#     c_je = y_intercept(m_je, (jx, jy))
#     mp_je = -1.0 / m_je
#     cp_je_d = y_intercept(mp_je, (dx, dy))
#     lx, ly = lines_intersect(m_je, c_je, mp_je, cp_je_d)
#     cp_je_k = y_intercept(mp_je, (kx, ky))
#     mx, my = lines_intersect(m_je, c_je, mp_je, cp_je_k)

 
#     green = [(lx,ly), (ex,ey), (xs[3],ys[3]), (dx,dy), (lx, ly)]
#     red = [(xs[1],ys[1]), (jx,jy), (lx,ly), (dx, dy), (xs[1], ys[1])]
#     yellow = [(kx,ky), (xs[2],ys[2]), (ex,ey), (mx, my), (kx,ky)]
#     blue = [(jx,jy), (kx,ky), (mx,my), (jx, jy)]
#     make_gif(side, red, green, blue, yellow)
# #     g_xs, g_ys = [p[1] for p in green], [p[2] for p in green]
# #     plot(xs, ys, legend=false, aspect_ratio=:equal)
# #     plot!(g_xs, g_ys, fill=(0, :lightgreen))
# #    # plot!([xs[1], jx, lx, dx, xs[1]], [ys[1], jy, ly, dy, ys[1]], fill=(0, :red))
# #     # plot!([jx, kx, mx, jx], [jy, ky, my, jy], fill=(0, :blue))
# #     # plot!([kx, xs[2], ex, mx, kx], [ky, ys[2], ey, my, ky], fill=(0, :yellow))
# #     red = rotate(red, (dx, dy), -π)
# #     r_xs, r_ys = [p[1] for p in red], [p[2] for p in red]
# #     plot!(r_xs, r_ys, fill=(0, :red))
# #     yv = rotate(yellow, 
# #     (ex, ey), π)
# #     y_xs, y_ys = [p[1] for p in yv], [p[2] for p in yv]
# #     plot!(y_xs, y_ys, fill=(0, :yellow))
# #     bv = rotate(blue, (ex, ey), π)
# #     b_xs, b_ys = [p[1] for p in bv], [p[2] for p in bv]
# #     # plot!(b_xs, b_ys, fill=(0, :blue))
# #     bv2 = rotate(zip(b_xs, b_ys), (b_xs[2], b_ys[2]), π)
# #     b2_xs, b2_ys = [p[1] for p in bv2], [p[2] for p in bv2]
# #     plot!(b2_xs, b2_ys, fill=(0, :blue))
# #     plot!([dx], [dy], marker=:circle)
# #     plot!([ex], [ey], marker=:circle)
# #     plot!([kx], [ky], marker=:circle)
# #     # plot!([fx], [fy], marker=:circle)
# #     # plot!([gx], [gy], marker=:circle)
# #     # plot!([hx], [hy], marker=:circle)
# #     # plot!([jx], [jy], marker=:circle)
# #     # plot!(xs_ae, ys_ae)
# #     # plot!(xs_bf, ys_bf)
# #     # plot!(xs_af, ys_af)
# #     # plot!(xs_hj, ys_hj)
# #     # plot!([xs[3], hx], [ys[3], hy])
# #     # plot!([jx, ex], [jy, ey])
# #     # plot!([dx, lx], [dy, ly])
# #     # plot!([kx, mx], [ky, my])
# end

end # module
