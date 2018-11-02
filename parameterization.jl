module ParamTest

import Base: +, -, *, ^, zero, one
import Graphics
using Gtk

sides = 5
delta = 0.0
resolution = 15
density = 0.1

# 3 => 0.0
# 4 => 0.0
# 5 => 0.4
# 6 => 0.35
# 7 => 0.2?

curves = []

struct Poly2D
    n::Int
    coeff::Array{Float64,2}
end

function +(u::Poly2D, v::Poly2D)
    if u.n < v.n
        m = copy(v.coeff)
        m[1:u.n,1:u.n] += u.coeff
        return Poly2D(v.n, m)
    elseif u.n > v.n
        m = copy(u.coeff)
        m[1:v.n,1:v.n] += v.coeff
        return Poly2D(u.n, m)
    end
    Poly2D(u.n, u.coeff + v.coeff)
end

-(u::Poly2D) = Poly2D(u.n, -copy(u.coeff))

-(u::Poly2D, v::Poly2D) = u + (-v)

function simplify(u::Poly2D)
    n = u.n
    for i in 1:n
        (u.coeff[i,end] != 0 || u.coeff[end,i] != 0) && return u
    end
    n == 1 ? u : simplify(Poly2D(n - 1, u.coeff[1:n-1,1:n-1]))
end

function *(u::Poly2D, v::Poly2D)
    n = u.n + v.n - 1
    m = zeros(Float64, n, n)
    len = u.n - 1
    for i in 1:v.n, j in 1:v.n
        x = v.coeff[i,j]
        if x != 0
            m[i:i+len,j:j+len] += u.coeff * x
        end
    end
    simplify(Poly2D(n, m))
end

*(u::Poly2D, x::Float64) = Poly2D(u.n, copy(u.coeff) * x)

function ^(u::Poly2D, n::Integer)
    @assert n >= 0 "Only non-negative (integer) exponents are supported."
    n == 0 && return one(Poly2D)
    n == 1 ? u : u ^ (n-1) * u
end

zero(::Type{Poly2D}) = Poly2D(1, zeros(Float64, 1, 1))

one(::Type{Poly2D}) = Poly2D(1, ones(Float64, 1, 1))

regularpoly(n) = [[0.5+cos(a)/2, 0.5+sin(a)/2] for a in range(0.0, length=n+1, stop=2pi)][1:n]

function line(p, q)
    m = zeros(Float64, 2, 2)
    d = q - p
    if abs(d[2]) > abs(d[1])
        x = d[1] / d[2]
        m[2,1] = 1
        m[1,2] = -x
        m[1,1] = -p[1] + x * p[2]
    else
        x = d[2] / d[1]
        m[2,1] = -x
        m[1,2] = 1
        m[1,1] = x * p[1] - p[2]
    end
    if m[2,1] + m[1,2] + 2 * m[1,1] < 0
        m *= -1
    end
    Poly2D(2, m)
end

lines(poly) = [line(p, q) for (p, q) in zip([[poly[end]]; poly[1:end-1]], poly)]

function polyeval(poly, uv)
    result = 0
    for i in 1:poly.n, j in 1:poly.n
        result += poly.coeff[i,j] * uv[1] ^ (i - 1) * uv[2] ^ (j - 1)
    end
    result
end

function normalized_lines(poly)
    n = length(poly)
    result = []
    for i in 1:n
        l = line(poly[mod1(i-1,n)], poly[i])
        l *= 1.0 / polyeval(l, poly[mod1(i+1,n)])
        push!(result, l)
    end
    result
end

affine_combine(p, x, q) = p * (1 - x) + q * x

function vertices(poly, resolution)
    n = length(poly)
    lines = [(poly[mod1(i-1,n)], poly[i]) for i in 1:n]
    center = [0.5, 0.5]
    result = [center]
    for j in 1:resolution
        coeff = j / resolution
        for k in 1:n, i in 0:j-1
            lp = affine_combine(lines[k][1], i / j, lines[k][2])
            push!(result, affine_combine(center, coeff, lp))
        end
    end
    result
end

function triangles(n, resolution)
    result = []
    inner_start = 1
    outer_vert = 2
    for layer in 1:resolution
        inner_vert = inner_start
        outer_start = outer_vert
        for side in 1:n
            vert = 1
            while true
                next_vert = side == n && vert == layer ? outer_start : outer_vert + 1
                push!(result, [inner_vert, outer_vert, next_vert])
                outer_vert += 1
                vert += 1
                vert == layer + 1 && break
                inner_next = side == n && vert == layer ? inner_start : inner_vert + 1
                push!(result, [inner_vert, next_vert, inner_next])
                inner_vert = inner_next
            end
        end
        inner_start = outer_start
    end
    result
end

function generate_curves()
    n = sides
    poly = regularpoly(n)
    L = normalized_lines(poly)
    H = []
    for i in 1:n
        im = mod1(i - 1, n)
        ip = mod1(i + 1, n)
        push!(H, L[i] * (one(Poly2D) - L[im] * L[ip] * delta))
    end
    vert = vertices(poly, resolution)
    tri = triangles(n, resolution)
    f(uv) = polyeval(H[1], uv)
    function slice(i, j)
        x, y = f(vert[i]), f(vert[j])
        if y > x
            x, y = y, x
            i, j = j, i
        end
        q1 = div(x, density)
        q2 = div(y, density)
        q1 - q2 == 1 && affine_combine(vert[j], (density * q1 - y) / (x - y), vert[i])
    end
    global curves = []
    foreach(tri) do t
        ab = slice(t[1], t[2])
        ac = slice(t[1], t[3])
        bc = slice(t[2], t[3])
        lst = filter(x -> x != false, [ab, ac, bc])
        if length(lst) == 2
            push!(curves, (lst[1], lst[2]))
        end
    end
end

function draw_polygon(ctx, poly, closep = false)
    if isempty(poly)
        return
    end
    Graphics.new_path(ctx)
    Graphics.move_to(ctx, poly[1][1], poly[1][2])
    for p in poly[2:end]
        Graphics.line_to(ctx, p[1], p[2])
    end
    if closep && length(poly) > 2
        Graphics.line_to(ctx, poly[1][1], poly[1][2])
    end
    Graphics.stroke(ctx)
end

function draw_segments(ctx, segments)
    for s in segments
        Graphics.new_path(ctx)
        Graphics.move_to(ctx, s[1][1], s[1][2])
        Graphics.line_to(ctx, s[2][1], s[2][2])
        Graphics.stroke(ctx)
    end
end

@guarded function draw_callback(canvas)
    ctx = Graphics.getgc(canvas)
    width = Graphics.width(canvas)
    height = Graphics.height(canvas)
    size = min(width, height)

    # White background
    Graphics.rectangle(ctx, 0, 0, width, height)
    Graphics.set_source_rgb(ctx, 1, 1, 1)
    Graphics.fill(ctx)

    # Input polygon
    Graphics.set_source_rgb(ctx, 1, 0, 0)
    Graphics.set_line_width(ctx, 2.0)
    draw_polygon(ctx, map(p -> p * size, regularpoly(sides)), true)

    # Generated curves
    Graphics.set_source_rgb(ctx, 0, 0, 0)
    Graphics.set_line_width(ctx, 1.0)
    draw_segments(ctx, map(s -> (s[1] * size, s[2] * size), curves))
end

function run()
    win = GtkWindow("Parameterization Test")
    vbox = GtkBox(:v)

    canvas = GtkCanvas(500, 500)
    draw(draw_callback, canvas)

    push!(win, vbox)
    push!(vbox, canvas)
    hbox = GtkBox(:h)
    push!(vbox, hbox)

    push!(hbox, GtkLabel("# of sides: "))
    sides_sb = GtkSpinButton(3, 10, 1)
    set_gtk_property!(sides_sb, :value, sides)
    signal_connect(sides_sb, "value-changed") do sb
        global sides = get_gtk_property(sb, :value, Int)
        generate_curves()
        draw(canvas)
    end
    push!(hbox, sides_sb)

    push!(hbox, GtkLabel("Delta: "))
    delta_sb = GtkSpinButton(0.0, 2.0, 0.01)
    set_gtk_property!(delta_sb, :value, delta)
    signal_connect(delta_sb, "value-changed") do sb
        global delta = get_gtk_property(sb, :value, Float64)
        generate_curves()
        draw(canvas)
    end
    push!(hbox, delta_sb)

    push!(hbox, GtkLabel("Resolution: "))
    resolution_sb = GtkSpinButton(5, 100, 5)
    set_gtk_property!(resolution_sb, :value, resolution)
    signal_connect(resolution_sb, "value-changed") do sb
        global resolution = get_gtk_property(sb, :value, Int)
        generate_curves()
        draw(canvas)
    end
    push!(hbox, resolution_sb)

    push!(hbox, GtkLabel("Density: "))
    density_sb = GtkSpinButton(0.01, 0.2, 0.05)
    set_gtk_property!(density_sb, :value, density)
    signal_connect(density_sb, "value-changed") do sb
        global density = get_gtk_property(sb, :value, Float64)
        generate_curves()
        draw(canvas)
    end
    push!(hbox, density_sb)

    showall(win)
    generate_curves()
    draw(canvas)
end

end # module
