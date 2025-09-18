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

function grooved_ellipsoid(p, a, b, c, d, k=50)
    smooth_max(ellipsoid(p, a, b, c), -shifted_ellipsoid(p, a, b, c, d), k)
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

function gen_mesh(body::EllipsoidalGroovedBody)
    @unpack a, b, c, groove_center = body
    mesh_from_function(
        p -> grooved_ellipsoid(p, a, b, c, groove_center), 
        origin=Vec3f(-a, -b, -c), 
        widths=Vec3f(2a, 2b, 2c),   
        samples=(150,150,150)
    )
end 