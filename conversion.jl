
module TensorKato

import Base: +, -, *, ^, zero, one, getindex, setindex!
using LinearAlgebra

exponent = 2
delta = 0.5

const Real = Float64 # BigFloat

struct Poly2D
    n::Int
    coeff::Array{Real,2}
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
    m = zeros(Real, n, n)
    len = u.n - 1
    for i in 1:v.n, j in 1:v.n
        x = v.coeff[i,j]
        if x != 0
            m[i:i+len,j:j+len] += u.coeff * x
        end
    end
    simplify(Poly2D(n, m))
end

*(u::Poly2D, x::Real) = Poly2D(u.n, copy(u.coeff) * x)
if Real != Float64
    *(u::Poly2D, x::Float64) = u * Real(x)
end

function ^(u::Poly2D, n::Integer)
    @assert n >= 0 "Only non-negative (integer) exponents are supported."
    n == 0 && return one(Poly2D)
    n == 1 ? u : u ^ (n-1) * u
end

zero(::Type{Poly2D}) = Poly2D(1, zeros(Real, 1, 1))

one(::Type{Poly2D}) = Poly2D(1, ones(Real, 1, 1))

"""
    regularpoly(n)

A regular `n`-gon on the inscribed circle of [0,1]x[0,1], consisting of an array of points.
For `n = 4` it returns the whole [0,1]x[0,1] square.
"""
function regularpoly(n)
    if n == 4
        [[0., 0], [1, 0], [1, 1], [0, 1]]
    else
        [[0.5+cos(a)/2, 0.5+sin(a)/2] for a in range(0.0, length=n+1, stop=2pi)][1:n]
    end
end

function line(p, q)
    m = zeros(Real, 2, 2)
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

binom(n, k) = Real(binomial(Int64(n), Int64(k)))

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
                    term *= L[im] ^ j * L[ip] ^ (ds - j) * binom(ds, j)
                    for l in 1:n
                        l == i && continue
                        lm = mod1(l - 1, n)
                        lp = mod1(l + 1, n)
                        term *= (L[lm] + L[lp]) ^ ds
                    end
                    term *= L[i] ^ k * (one(Poly2D) - L[i]) ^ (dh - k) * binom(dh, k)
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
    n == 3 && return dgregory(ribbons)
    poly = regularpoly(n)
    L = normalized_lines(poly)
    function multiplier(on_0, on_d1)
        term = one(Poly2D)
        for i in 1:n
            i in on_0 && continue
            im = mod1(i - 1, n)
            ip = mod1(i + 1, n)
            d = i in on_d1 ? ds - 1 : ds
            term *= (L[im] + L[ip]) ^ d
        end
        term
    end
    denominator = sum(i -> blend2(L, i), 1:n) * multiplier([], [])
    numerator = [zero(Poly2D), zero(Poly2D), zero(Poly2D)]
    for c in 1:3
        for i in 1:n
            im = mod1(i - 1, n)
            imm = mod1(im - 1, n)
            ip = mod1(i + 1, n)
            function ribbon(i, i1, s, s1, h, h1)
                R = zero(Poly2D)
                for j in 0:ds
                    p1 = ribbons.cpts[i-1,j,0][c]
                    p2 = ribbons.cpts[i-1,j,1][c]
                    term = (h + h1) * p1 + h * (p2 - p1) * float(ds)
                    term *= s ^ j * s1 ^ (ds - j) * binom(ds, j)
                    R += term
                end
                R * multiplier([i], [i1])
            end
            r1 = ribbon(im, i, L[imm], L[i], L[im], L[ip])
            r2 = ribbon(i, im, L[im], L[ip], L[i], L[imm])
            corner = ribbons.cpts[i-1,0,0][c]
            twist  = ribbons.cpts[i-1,1,1][c]
            left   = ribbons.cpts[i-1,1,0][c]
            right  = ribbons.cpts[i-1,0,1][c]
            q = multiplier([], []) * corner +
                L[im] * multiplier([], [i]) * float(ds) * (left - corner) +
                L[i] * multiplier([], [im]) * float(ds) * (right - corner) +
                L[im] * L[i] * multiplier([], [i, im]) * float(ds * ds) *
                (twist - left - right + corner)
            numerator[c] += (r1 + r2 - q) * blend2(L, i)
        end
    end
    (numerator, denominator)
end

function dgregory(ribbons)
    n, ds = ribbons.n, ribbons.d
    poly = regularpoly(n)
    L = normalized_lines(poly)
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
                    term *= s ^ j * (one(Poly2D) - s) ^ (ds - j) * binom(ds, j)
                    R += term
                end
                R
            end
            r1 = ribbon(im, one(Poly2D) - L[i], L[im])
            r2 = ribbon(i, L[im], L[i])
            corner = ribbons.cpts[i-1,0,0][c]
            twist  = ribbons.cpts[i-1,1,1][c]
            left   = ribbons.cpts[i-1,1,0][c]
            right  = ribbons.cpts[i-1,0,1][c]
            q = one(Poly2D) * corner +
                L[im] * float(ds) * (left - corner) +
                L[i] * float(ds) * (right - corner) +
                L[im] * L[i] * float(ds * ds) * (twist - left - right + corner)
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
    result = zeros(Real, d + 1, d + 1)
    for i in 0:d
        for j in 0:d-i
            result[i+1,i+j+1] = binom(d, i) * binom(d - i, j) * (isodd(j) ? -1 : 1)
        end
    end
    result
end

function tobezier(d)            # = frombezier(d) ^ -1
    result = zeros(Real, d + 1, d + 1)
    for i in 0:d
        for j in i:d
            result[i+1,j+1] = binom(j, i)
        end
        result[i+1,:] /= binom(d, i)
    end
    result
end

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
    bernstein(k, x) = binom(d, k) * x ^ k * (1 - x) ^ (d - k)
    result = [0, 0, 0, 0]
    for i in 0:d, j in 0:d
        result += surf[i+1,j+1,:] * bernstein(i, uv[1]) * bernstein(j, uv[2])
    end
    result[1:3] / result[4]
end

function bezier(ribbons)
    d = ribbons.d
    result = Array{Float64}(undef, d + 1, d + 1, 4)
    center = d / 2
    for i in 0:d, j in 0:d
        if i < center
            result[i+1,j+1,:] = [ribbons.cpts[3,d-j,i]; 1]
        elseif i > center
            result[i+1,j+1,:] = [ribbons.cpts[1,j,d-i]; 1]
        elseif j < center
            result[i+1,j+1,:] = [ribbons.cpts[0,i,j]; 1]
        elseif j > center
            result[i+1,j+1,:] = [ribbons.cpts[2,d-i,d-j]; 1]
        else # i == j == center
            result[i+1,j+1,:] = [ribbons.center; 1]
        end
    end
    result
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
    center :: Point
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
        l = Int(floor((d + 1) / 2))
        cp = 1 + Int(floor(d / 2))
        cp = n * cp * l + 1
        side, col, row = 0, 0, 0
        center = read_numbers(f, Float64)
        result = BezierPatch(n, d, center, Dict())
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
    bernstein(n, k, x) = binomial(n, k) * x ^ k * (1 - x) ^ (n - k)
    samples = [[u, v] for u in range(0.0, stop=1.0, length=resolution)
                      for v in range(0.0, stop=1.0, length=resolution)]
    d = ribbons.d
    verts = map(samples) do p
        result = [0, 0, 0]
        for i in 0:d, j in 0:1
            result += ribbons[index,i,j] * bernstein(d, i, p[1]) * bernstein(1, j, p[2])
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
    if n == 4
        surf = bezier(ribbons)
        write_cnet(surf, "$filename-cnet.obj")
        write_surface(evaltensor, surf, 0, resolution, "$filename.obj")
        return
    end
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


# S-patch

const Index = Vector{Int}

struct SPatch
    n :: Int
    d :: Int
    cpts :: Dict{Index,Point}
end

getindex(s::SPatch, si::Index) = s.cpts[si]
setindex!(s::SPatch, p::Point, si::Index) = s.cpts[si] = p
get!(s::SPatch, si::Index, p::Point) = get!(s.cpts, si, p)

indices(n, d) = n == 1 ? [d] : [[i; si] for i in 0:d for si in indices(n - 1, d - i)]

multinomial(si) = factorial(sum(si)) ÷ prod(map(factorial, si))

multibernstein(si, bc) = prod(map(^, bc, si)) * Real(multinomial(si))

function read_spatch(filename)
    read_numbers(lst, numtype) = map(s -> parse(numtype, s), lst)
    open(filename) do f
        n, d = read_numbers(split(readline(f)), Int)
        size = binom(n + d - 1, d)
        result = SPatch(n, d, Dict())
        for _ in 1:size
            line = split(readline(f))
            index = read_numbers(line[1:n], Int)
            point = read_numbers(line[n+1:end], Float64)
            result.cpts[index] = point
        end
        result
    end
end

function spatch(surf)
    n, d = surf.n, surf.d
    poly = regularpoly(n)
    L = normalized_lines(poly)
    bc = [prod(j -> j == i || j == mod1(i - 1, n) ? one(Poly2D) : L[j], 1:n) for i in 1:n]
    denominator = n <= 4 ? one(Poly2D) : sum(bc) ^ d
    if length(iterate(surf.cpts)[1][2]) == 4
        denominator *= sum(p -> multibernstein(p[1], bc) * p[2][4], collect(surf.cpts))
    end
    numerator = [zero(Poly2D), zero(Poly2D), zero(Poly2D)]
    for (i, q) in surf.cpts, c in 1:3
        numerator[c] += multibernstein(i, bc) * q[c]
    end
    (numerator, denominator)
end

function spatch_test(filename, resolution = 15)
    surf = read_spatch("$filename.sp")
    tsurf = tensor(spatch(surf))
    write_cnet(tsurf, "$filename-spatch-tensor-cnet.obj")
    write_surface(evaltensor, tsurf, surf.n, resolution, "$filename-spatch-tensor.obj")
    write_surface(evaltensor, tsurf, 0, resolution, "$filename-spatch-tensor-full.obj")
end


# Warren's patch

function warren(surf)
    n, d = surf.n, surf.d
    @assert n == 3 "Warren's patch works only on Bezier triangles"
    u = Poly2D(2, [0. 1; 0 -1])
    v = Poly2D(2, [0. 0; 0 1])
    w = Poly2D(2, [1. -1; 0 0])
    bc = [u, v, w]
    denominator = one(Poly2D)
    if length(iterate(surf.cpts)[1][2]) == 4
        denominator *= sum(p -> multibernstein(p[1], bc) * p[2][4], collect(surf.cpts))
    end
    numerator = [zero(Poly2D), zero(Poly2D), zero(Poly2D)]
    for (i, q) in surf.cpts, c in 1:3
        numerator[c] += multibernstein(i, bc) * q[c]
    end
    (numerator, denominator)
end

function warren_test(filename, resolution = 15)
    surf = read_spatch("$filename.sp")
    tsurf = tensor(warren(surf))
    write_cnet(tsurf, "$filename-warren-tensor-cnet.obj")
    write_surface(evaltensor, tsurf, 0, resolution, "$filename-warren-tensor.obj")
end

end # module
