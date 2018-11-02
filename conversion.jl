module TensorKato

import Base: +, -, *, ^, zero, one, getindex, setindex!
using LinearAlgebra

exponent = 2
delta = 0.5

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

blend(L, i) = prod(j -> L[j] ^ exponent, setdiff(1:length(L), i))

binom(n, k) = binomial(Int64(n), Int64(k))

function kato(ribbons)
    n, ds = ribbons.n, ribbons.d
    dh = 1
    poly = regularpoly(n)
    L = lines(poly)
    denominator = sum(i -> blend(L, i), 1:n)
    for i in 1:n
        im = mod1(i - 1, n)
        ip = mod1(i + 1, n)
        denominator *= (L[im] + L[ip]) ^ ds
    end
    numerator = [zero(Poly2D), zero(Poly2D), zero(Poly2D)]
    for c in 1:3
        for i in 1:n
            im = mod1(i - 1, n)
            ip = mod1(i + 1, n)
            R = zero(Poly2D)
            for j in 0:ds
                for k in 0:dh
                    p = ribbons.cpts[i-1,j,k][c]
                    if k == 1
                        # Scale ribbons
                        p0 = ribbons.cpts[i-1,j,0][c]
                        p = p0 + (p - p0) * ds
                    end
                    term = one(Poly2D) * p
                    term *= L[im] ^ j * L[ip] ^ (ds - j) * float(binom(ds, j))
                    for l in 1:n
                        l == i && continue
                        lm = mod1(l - 1, n)
                        lp = mod1(l + 1, n)
                        term *= (L[lm] + L[lp]) ^ ds
                    end
                    term *= L[i] ^ k * (one(Poly2D) - L[i]) ^ (dh - k) * float(binom(dh, k))
                    R += term
                end
            end
            numerator[c] += R * blend(L, i)
        end
    end
    (numerator, denominator)
end

blend2(L, i) = prod(j -> L[j] ^ 2, setdiff(1:length(L), [mod1(i-1,length(L)), i]))

function gregory(ribbons)
    n, ds = ribbons.n, ribbons.d
    poly = regularpoly(n)
    L = normalized_lines(poly)
    H = []
    for i in 1:n
        im = mod1(i - 1, n)
        ip = mod1(i + 1, n)
        push!(H, L[i] * (one(Poly2D) - L[im] * L[ip] * delta))
    end
    denominator = sum(i -> blend2(L, i), 1:n)
    numerator = [zero(Poly2D), zero(Poly2D), zero(Poly2D)]
    for c in 1:3
        for i in 1:n
            im = mod1(i - 1, n)
            function ribbon(i, s, h)
                R = zero(Poly2D)
                for j in 0:ds
                    p1 = ribbons.cpts[i-1,j,0][c]
                    p2 = ribbons.cpts[i-1,j,1][c]
                    term = one(Poly2D) * p1 + h * (p2 - p1) * float(ds)
                    term *= s ^ j * (one(Poly2D) - s) ^ (ds - j) * float(binom(ds, j))
                    R += term
                end
                R
            end
            r1 = ribbon(im, one(Poly2D) - H[i], H[im])
            r2 = ribbon(i, H[im], H[i])
            corner = ribbons.cpts[i-1,0,0][c]
            twist  = ribbons.cpts[i-1,1,1][c]
            left   = ribbons.cpts[i-1,1,0][c]
            right  = ribbons.cpts[i-1,0,1][c]
            q = one(Poly2D) * corner +
                H[im] * float(ds) * (left - corner) +
                H[i] * float(ds) * (right - corner) +
                H[im] * H[i] * float(ds * ds) * (twist - left - right + corner)
            numerator[c] += (r1 + r2 - q) * blend2(L, i)
        end
    end
    (numerator, denominator)
end

function eval(kato, uv)
    num, den = kato
    d = 0
    for i in 1:den.n
        u = uv[1] ^ (i - 1)
        for j in 1:den.n
            v = uv[2] ^ (j - 1)
            d += den.coeff[i,j] * u * v
        end
    end
    p = [0., 0, 0]
    for c in 1:3
        for i in 1:num[c].n
            u = uv[1] ^ (i - 1)
            for j in 1:num[c].n
                v = uv[2] ^ (j - 1)
                p[c] += num[c].coeff[i,j] * u * v
            end
        end
    end
    p / d
end

function frombezier(d)
    result = zeros(d + 1, d + 1)
    for i in 0:d
        for j in 0:d-i
            result[i+1,i+j+1] = binom(d, i) * binom(d - i, j) * (isodd(j) ? -1 : 1)
        end
    end
    result
end

tobezier(d) = frombezier(d) ^ -1

function tobezier(coeff, d)
    C = tobezier(d)
    d1 = size(coeff, 1)
    M = coeff
    if d1 < d + 1
        M = zeros(d + 1, d + 1)
        M[1:d1,1:d1] = coeff
    end
    C' * M * C
end

function tensor(kato)
    num, den = kato
    println("Degree: $(num[1].n-1)/$(den.n-1)")
    d = max(num[1].n, den.n) - 1
    w = tobezier(den.coeff, d)
    px = tobezier(num[1].coeff, d)
    py = tobezier(num[2].coeff, d)
    pz = tobezier(num[3].coeff, d)
    result = Array{Float64}(undef, d + 1, d + 1, 4)
    for i in 0:d, j in 0:d
        result[i+1,j+1,1] = px[i+1,j+1]
        result[i+1,j+1,2] = py[i+1,j+1]
        result[i+1,j+1,3] = pz[i+1,j+1]
        result[i+1,j+1,4] = w[i+1,j+1]
    end
    result
end

function evaltensor(surf, uv)
    d = size(surf, 1) - 1
    bezier(k, x) = binom(d, k) * x ^ k * (1 - x) ^ (d - k)
    result = [0, 0, 0, 0]
    for i in 0:d, j in 0:d
        result += surf[i+1,j+1,:] * bezier(i, uv[1]) * bezier(j, uv[2])
    end
    result[1:3] / result[4]
end


# I/O

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

const Point = Vector{Float64}
const GBIndex = NTuple{3, Int}

struct BezierPatch
    n :: Int
    d :: Int
    cpts :: Dict{GBIndex,Point}
end

getindex(s::BezierPatch, idx::GBIndex) = s.cpts[idx]
getindex(s::BezierPatch, i::Int, j::Int, k::Int) = s.cpts[i,j,k]
setindex!(s::BezierPatch, v::Point, idx::GBIndex) = s.cpts[idx] = v
setindex!(s::BezierPatch, v::Point, i::Int, j::Int, k::Int) = s.cpts[i,j,k] = v

function read_ribbons(filename)
    read_numbers(f, numtype) = map(s -> parse(numtype, s), split(readline(f)))
    local result
    open(filename) do f
        n, d = read_numbers(f, Int)
        result = BezierPatch(n, d, Dict())
        l = Int(floor((d + 1) / 2))
        cp = 1 + Int(floor(d / 2))
        cp = n * cp * l + 1
        side, col, row = 0, 0, 0
        readline(f)
        for i in 1:cp-1
            if col >= d - row
                side += 1
                if side >= n
                    side = 0
                    row += 1
                end
                col = row
            end
            p = read_numbers(f, Float64)
            result[side,col,row] = p
            if col < l
                result[mod(side-1,n),d-row,col] = p
            elseif d - col < l
                result[mod(side+1,n),row,d-col] = p
            end
            col += 1
        end
    end
    result
end

function writeOBJ(verts, tris, filename)
    open(filename, "w") do f
        for v in verts
            println(f, "v $(v[1]) $(v[2]) $(v[3])")
        end
        for t in tris
            println(f, "f $(t[1]) $(t[2]) $(t[3])")
        end
    end
end

function write_frames(ribbons, filename)
    g1 = filter(idx -> idx[3] <= 1, keys(ribbons.cpts))
    mapping = Dict(map(p -> (p[2], p[1]), enumerate(g1)))
    open(filename, "w") do f
        for idx in g1
            p = ribbons[idx]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
        end
        for idx in g1
            i, j, k = idx
            from = mapping[idx]
            for next in [(i, j + 1, k), (i, j, k + 1)]
                !haskey(mapping, next) && continue
                to = mapping[next]
                println(f, "l $from $to")
            end
        end
    end
end

function write_cnet(surf, filename)
    d = size(surf, 1) - 1
    open(filename, "w") do f
        for i in 0:d, j in 0:d
            p = surf[i+1,j+1,1:3] / surf[i+1,j+1,4]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
        end
        for i in 1:d+1, j in 1:d
            index = (i - 1) * (d + 1) + j
            println(f, "l $index $(index+1)")
            index = (j - 1) * (d + 1) + i
            println(f, "l $index $(index+d+1)")
        end
    end
end

# - fn is the evaluator function, n the # of sides
# - uses the whole [0,1]x[0,1] domain when n = 0
function write_surface(fn, surf, n, resolution, filename)
    local verts, tris
    if n == 0
        verts = [fn(surf, [u, v])
                 for u in range(0.0, stop=1.0, length=resolution)
                 for v in range(0.0, stop=1.0, length=resolution)]
        tris = []
        for i in 2:resolution, j in 2:resolution
            index = (j - 1) * resolution + i
            push!(tris, [index - resolution - 1, index - resolution, index])
            push!(tris, [index, index - 1, index - resolution - 1])
        end
    else
        poly = regularpoly(n)
        verts = [fn(surf, p) for p in vertices(poly, resolution)]
        tris = triangles(n, resolution)
    end
    writeOBJ(verts, tris, filename)
end

function write_ribbon(ribbons, filename, index, resolution)
    bezier(n, k, x) = binomial(n, k) * x ^ k * (1 - x) ^ (n - k)
    samples = [[u, v] for u in range(0.0, stop=1.0, length=resolution)
                      for v in range(0.0, stop=1.0, length=resolution)]
    d = ribbons.d
    verts = map(samples) do p
        result = [0, 0, 0]
        for i in 0:d, j in 0:1
            result += ribbons[index,i,j] * bezier(d, i, p[1]) * bezier(1, j, p[2])
        end
        result
    end
    tris = []
    for i in 2:resolution, j in 2:resolution
        index = (j - 1) * resolution + i
        push!(tris, [index - resolution - 1, index - resolution, index])
        push!(tris, [index, index - 1, index - resolution - 1])
    end
    writeOBJ(verts, tris, filename)
end


# Test

function test(filename, resolution = 15, surftype = :gregory)
    ribbons = read_ribbons("$filename.gbp")
    n = ribbons.n
    surf = surftype == :gregory ? gregory(ribbons) : kato(ribbons)
    tsurf = tensor(surf)
    write_frames(ribbons, "$filename-frames.obj")
    write_surface(eval, surf, n, resolution, "$filename.obj")
    write_surface(eval, surf, 0, resolution, "$filename-full.obj")
    write_cnet(tsurf, "$filename-cnet.obj")
    write_surface(evaltensor, tsurf, n, resolution, "$filename-tensor.obj")
    write_surface(evaltensor, tsurf, 0, resolution, "$filename-tensor-full.obj")
    for i in 1:n
        write_ribbon(ribbons, "$filename-$i.obj", i - 1, resolution)
    end
end

end # module
