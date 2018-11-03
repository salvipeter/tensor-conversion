module ParamTest

import Graphics
using Gtk

using ..TensorKato: Poly2D, regularpoly, normalized_lines,
    vertices, triangles, polyeval, affine_combine

sides = 5
delta = 0.0
resolution = 15
density = 0.1

curves = []

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
    # f(uv) = polyeval(H[1], uv)
    f(uv) = polyeval(L[1], uv) / polyeval(L[1] + L[3], uv)
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

macro add_spinbox(layout, label, var, min, max, step, valtype)
    return quote
        push!($(esc(layout)), GtkLabel($label))
        local sb = GtkSpinButton($min, $max, $step)
        set_gtk_property!(sb, :value, $var)
        signal_connect(sb, "value-changed") do sb
            global $var = get_gtk_property(sb, :value, $valtype)
            generate_curves()
            draw($(esc(:(canvas))))
        end
        push!($(esc(layout)), sb)
    end
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

    @add_spinbox(hbox, "# of sides: ", sides, 3, 10, 1, Int)
    @add_spinbox(hbox, "Delta: ", delta, 0.0, 2.0, 0.01, Float64)
    @add_spinbox(hbox, "Resolution: ", resolution, 5, 100, 5, Int)
    @add_spinbox(hbox, "Density: ", density, 0.01, 0.2, 0.05, Float64)

    showall(win)
    generate_curves()
    draw(canvas)
end

end # module
