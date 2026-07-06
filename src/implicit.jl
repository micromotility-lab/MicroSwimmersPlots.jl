function smooth_max(a, b, k)
    m = max(a, b)
    m + log(exp(k*(a - m)) + exp(k*(b - m))) / k
end
# smooth_max(a, b, k) = log(exp(k*a) + exp(k*b)) / k

function ellipsoid(p, a, b, c)
    x, y, z = p
    x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
end

function shifted_ellipsoid(p, a, b, c, d)
    x, y, z = p .- d  # subtract the shift vector d = (dx, dy, dz)
    x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
end

function shifted_rotated_ellipsoid(p, a, b, c, d, R)
    x, y, z = R' * (p .- d)  
    x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
end

function cylinder(p, a, b, h, d)
    x, y, z = p .- d
    max(x^2/a^2 + y^2/b^2 - 1, -z-h, z - h)
end

plane(p, p0, n) = dot(p .- p0, n)

function ellipsoidal_grooved_ellipsoid(p, a, b, c, g_a, g_b, g_c, d, R, k=50)
    smooth_max(ellipsoid(p, a, b, c), -shifted_rotated_ellipsoid(p, g_a, g_b, g_c, d, R), k)
end

function sliced_ellipsoidal_grooved_ellipsoid(
    p, a, b, c, g_a, g_b, g_c, d, R, p0, n, k=50
)
    body_with_groove = ellipsoidal_grooved_ellipsoid(p, a, b, c, g_a, g_b, g_c, d, R, k)

    max(body_with_groove, plane(p, p0, n))
end

function cylindrical_grooved_ellipsoid(p, a, b, c, g_a, g_b, h, d, k=50)
    smooth_max(ellipsoid(p, a, b, c), -cylinder(p, g_a, g_b, h, d), k)
end