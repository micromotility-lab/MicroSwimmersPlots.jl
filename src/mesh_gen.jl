using LinearAlgebra
using GeometryBasics


smooth_max(a, b, k) = log(exp(k*a) + exp(k*b)) / k

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


function mesh_from_function(f::Function; origin=Vec3f(0., 0., 0), widths=Vec3f(1., 1., 1.), samples=(150,150,150))
    corner = origin + widths
    xs = range(origin[1], corner[1], samples[1])
    ys = range(origin[2], corner[2], samples[2])
    zs = range(origin[3], corner[3], samples[3])
    points, faces = isosurface(
        [f(Vec3f(x,y,z)) for x in xs, y in ys, z in zs],
        MarchingCubes(),
        xs, ys, zs
    )
    GeometryBasics.Mesh(Point3f.(points), TriangleFace.(faces))
end


function gen_mesh(body::EllipsoidBody)
    @unpack a, b, c = body
    mesh_from_function(
        p -> sum(x -> x^2, p ./ [a, b, c]) - 1, 
        origin=Vec3f(-a, -b, -c),
        widths=Vec3f(2a, 2b, 2c),
        samples=(150,150,150)
    )
end 

function gen_mesh_sliced(body::EllipsoidalGroovedBody; p0=Point3f(0.), n=Vec3f(0., 0., 1.), samples=(150,150,150))
    @unpack a, b, c, g_a, g_b, g_c, groove_center, orientation = body

    mesh_from_function(
        p -> sliced_ellipsoidal_grooved_ellipsoid(
            p, a, b, c, g_a, g_b, g_c, groove_center, orientation, p0, n
        ),
        origin=Vec3f(-a, -b, -c),
        widths=Vec3f(2a, 2b, 2c),
        samples=samples
    )
end

function gen_mesh(body::EllipsoidalGroovedBody)
    @unpack a, b, c, g_a, g_b, g_c, groove_center, orientation = body
    mesh_from_function(
        p -> ellipsoidal_grooved_ellipsoid(p, a, b, c, g_a, g_b, g_c, groove_center, orientation), 
        origin=Vec3f(-a, -b, -c), 
        widths=Vec3f(2a, 2b, 2c),   
        samples=(150,150,150)
    )
end 

function gen_mesh(body::FlatGroovedBody)
    @unpack a, b, c, g_a, g_b, groove_floor_center = body
    g_depth = groove_floor_center[3]
    cyl_height = c - g_depth
    cyl_center = groove_floor_center .+ [0., 0., 0.5cyl_height]
    mesh_from_function(
        p -> cylindrical_grooved_ellipsoid(p, a, b, c, g_a, g_b, 0.5cyl_height, cyl_center), 
        origin=Vec3f(-a, -b, -c), 
        widths=Vec3f(2a, 2b, 2c), 
        samples=(150,150,150)
    )
end
