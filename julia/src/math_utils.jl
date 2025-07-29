module MathUtils

export integrate_romberg, brentq, gauss_legendre

using Roots

"""
    integrate_romberg(f, a, b; tol=1e-8, maxit=25, minsteps=0, abs_tol=false)

Integrate function `f` from `a` to `b` using Romberg integration.
"""
function integrate_romberg(f, a, b; tol=1e-8, maxit=25, minsteps=0, abs_tol=false)
    h = 0.5 * (b - a)
    gprev = h * (f(a) + f(b))
    g = [gprev]
    nint = 1
    err = Inf
    i = 0
    while true
        i += 1
        if i > maxit || (i > 5 && abs(err) < tol && nint > minsteps)
            break
        end
        g0 = zero(gprev)
        for k in 1:nint
            g0 += f(a + (2k - 1) * h)
        end
        g0 = 0.5 * gprev + h * g0
        h *= 0.5
        nint *= 2
        jmax = min(i, 5)
        fourj = 1.0
        for j in 1:jmax
            fourj *= 4.0
            g1 = g0 + (g0 - g[j]) / (fourj - 1.0)
            g[j] = g0
            g0 = g1
        end
        err = abs(abs_tol ? gprev - g0 : (abs(g0) > tol ? 1.0 - gprev / g0 : gprev))
        gprev = g0
        if length(g) < jmax + 1
            push!(g, g0)
        else
            g[jmax + 1] = g0
        end
    end
    return gprev
end

"""
    brentq(f, a, b; tol=1e-8)

Find a root of `f` in [`a`,`b`] using the `Roots` package implementation
of the Brent method.
"""
brentq(f, a, b; tol=1e-8) = find_zero((f, a, b), Brent(); xtol=tol)

"""
    gauss_legendre(n)

Return nodes and weights for Gauss-Legendre quadrature with `n` points.
"""
function gauss_legendre(n::Int)
    x = zeros(n)
    w = zeros(n)
    m = (n + 1) ÷ 2
    for i in 1:m
        z = cos(pi * (i - 0.25) / (n + 0.5))
        z1 = 0.0
        while abs(z - z1) > 1e-15
            p1 = 1.0
            p2 = 0.0
            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2*j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z*p1 - p2) / (z^2 - 1)
            z1 = z
            z -= p1 / pp
        end
        x[i] = -z
        x[n + 1 - i] = z
        w[i] = 2 / ((1 - z^2) * pp^2)
        w[n + 1 - i] = w[i]
    end
    return x, w
end

end # module
