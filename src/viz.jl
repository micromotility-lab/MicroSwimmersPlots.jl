const I3f = Mat3f(I)

abstract type PlotBuffer end
abstract type FlagellumBuffer <: PlotBuffer end

function update_pts!(out_vec::Vector{Point3f}, in_vec::Vector{Point3f}, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    @inbounds for i in eachindex(in_vec)
        out_vec[i] = location + orientation*in_vec[i]
    end 
end

function update_pts!(out_vec::Vector{Point3f}, in_vec::Matrix{<:Number}, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    @inbounds for i in axes(in_vec, 2)
        out_vec[i] = location + orientation*Point3f(@view in_vec[:,i])
    end 
end

function update_buffer_observable!(B::Observable{P}, fb::FluidBoundary) where {P <: PlotBuffer}
    update_buffer!(B[], fb)
    notify(B)
end


function viz(fb::FluidBoundary; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    B = viz!(ax, fb; kwargs...)
    display(fig)
    B
end

function viz!(ax, buf::PlotBuffer; kwargs...)
    B = Observable(buf)
    viz!(ax, B; kwargs...)
    B
end

function viz!(ax, fb::FluidBoundary; kwargs...)
    buf = get_buffer(fb)
    update_buffer!(buf, fb)
    viz!(ax, buf; kwargs...)
end


# Cell bodies

mutable struct CellBodyBuffer <: PlotBuffer
    mesh::GeometryBasics.Mesh
    ref_pts::Vector{Point3f}
end

function get_buffer(body::CellBody)
    mesh = gen_mesh(body.model)
    ref_pts = copy(coordinates(mesh))
    CellBodyBuffer(mesh, ref_pts)
end

function update_buffer!(buf::CellBodyBuffer, body::CellBody, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    update_pts!(buf.mesh.position, buf.ref_pts, Point3f(body.points.location), Mat3f(body.points.orientation))
    update_pts!(buf.mesh.position, buf.mesh.position, location, orientation)
end

function update_mesh!(buf::CellBodyBuffer, body::CellBody, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    buf.mesh = gen_mesh(body.model)
    buf.ref_pts = copy(coordinates(buf.mesh))
    update_pts!(buf.mesh.position, buf.ref_pts, Point3f(body.points.location), Mat3f(body.points.orientation))
    update_pts!(buf.mesh.position, buf.mesh.position, location, orientation)
end

viz!(ax, B::Observable{CellBodyBuffer}; kwargs...) = mesh!(ax, @lift($B.mesh); kwargs...)

# Flagellum 

struct BareFlagellumBuffer <: FlagellumBuffer
    pts::Vector{Point3f}
end

get_buffer(f::BareFlagellum) = BareFlagellumBuffer(Vector{Point3f}(undef, f.points.N))

function update_buffer!(buf::BareFlagellumBuffer, f::BareFlagellum, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    update_pts!(buf.pts, f.points.force_pts, Point3f(f.points.location), Mat3f(f.points.orientation))
    update_pts!(buf.pts, buf.pts, location, orientation)
end

function viz!(ax, B::Observable{BareFlagellumBuffer}; linewidth=3, color=:forestgreen)
    lines!(ax, @lift($B.pts), linewidth=linewidth, color=color)
end

struct VanedFlagellumBuffer <: FlagellumBuffer
    flagellum_pts::Vector{Point3f}
    vane_pts::Vector{Point3f}
    vane_x::Matrix{Float32}
    vane_y::Matrix{Float32}
    vane_z::Matrix{Float32}
end

get_buffer(vf::VanedFlagellum) = VanedFlagellumBuffer(
    Vector{Point3f}(undef, vf.N_f),
    Vector{Point3f}(undef, vf.N_v*(vf.N_height+1)),
    Matrix{Float32}(undef, vf.N_v, vf.N_height + 1),
    Matrix{Float32}(undef, vf.N_v, vf.N_height + 1),
    Matrix{Float32}(undef, vf.N_v, vf.N_height + 1)
)

function update_buffer!(buf::VanedFlagellumBuffer, vf::VanedFlagellum, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    update_pts!(buf.flagellum_pts, vf.points.force_pts[:, 1:vf.N_f], Point3f(vf.points.location), Mat3f(vf.points.orientation))
    update_pts!(buf.flagellum_pts, buf.flagellum_pts, location, orientation)
    update_pts!(buf.vane_pts, get_vane_pts(vf), Point3f(vf.points.location), Mat3f(vf.points.orientation))
    update_pts!(buf.vane_pts, buf.vane_pts, location, orientation)
    for i in eachindex(buf.vane_pts)
        buf.vane_x[i] = buf.vane_pts[i][1]
        buf.vane_y[i] = buf.vane_pts[i][2]
        buf.vane_z[i] = buf.vane_pts[i][3]
    end
end

function viz!(ax, B::Observable{VanedFlagellumBuffer}; linewidth=3, color=:forestgreen)
    lines!(ax, @lift($B.flagellum_pts), linewidth=linewidth, color=color)
    surface!(ax,
        @lift($B.vane_x), @lift($B.vane_y), @lift($B.vane_z),
        colormap=:RdPu
    )
end

# Flagellate

struct FlagellateBuffer{F <: FlagellumBuffer} <: PlotBuffer
    body_buffer::CellBodyBuffer
    flagella_buffers::Vector{F}
end

function get_buffer(flg::Flagellate)
    body_buf = get_buffer(flg.body)
    flagella_bufs = [get_buffer(f) for f in flg.flagella]
    FlagellateBuffer(body_buf, flagella_bufs)
end

function update_buffer!(buf::FlagellateBuffer, flg::Flagellate)
    update_buffer!(buf.body_buffer, flg.body, Point3f(flg.points.location), Mat3f(flg.points.orientation))
    # update_pts!(buf.body_buffer.mesh.position, buf.body_buffer.ref_pts, Point3f(flg.points.location), Mat3f(flg.points.orientation))
    for (i,f) in enumerate(flg.flagella)
        fb = buf.flagella_buffers[i]
        update_buffer!(fb, f, Point3f(flg.points.location), Mat3f(flg.points.orientation))
    end
end
  
function viz!(ax, B::Observable{FlagellateBuffer{T}}; bodycolor=Makie.wong_colors()[1], linewidth=3, color=:forestgreen) where {T <: FlagellumBuffer}
    mesh!(ax, @lift($B.body_buffer.mesh), color=bodycolor)
    for i in eachindex(B[].flagella_buffers)
        viz!(ax, @lift($B.flagella_buffers[i]); linewidth=linewidth, color=color)
    end
end

viz(prob::SwimmingProblem) = viz(prob.microswimmer)

function viz(prob::SwimmingTrajectoryProblem)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    lines!(ax, prob.traj.x, color=:red)
    B = viz!(ax, prob.swimming_problem.microswimmer)
    s = Slider(fig[2,1], range=1:length(prob.traj.t), startvalue=1)
    on(s.value) do val
        move_boundary!(prob, val)
        update_buffer_observable!(B, prob.swimming_problem.microswimmer)
    end
    fig
end


# Velocity fields



@inline function _proj2(plane::Symbol, x)
    plane === :xy && return Point2f(x[1], x[2])
    plane === :xz && return Point2f(x[1], x[3])
    plane === :yz && return Point2f(x[2], x[3])
    throw(ArgumentError("plane must be :xy, :xz, or :yz"))
end

# @inline function _slice3(plane::Symbol, a::Real, b::Real, c::Real)
#     # given (a,b) in plane coordinates and fixed c, return 3D point
#     plane === :xy && return (a, b, c)
#     plane === :xz && return (a, c, b)
#     plane === :yz && return (c, a, b)
#     throw(ArgumentError("plane must be :xy, :xz, or :yz"))
# end

@inline function _labels(plane::Symbol)
    plane === :xy && return (L"x", L"y")
    plane === :xz && return (L"x", L"z")
    plane === :yz && return (L"y", L"z")
    throw(ArgumentError("plane must be :xy, :xz, or :yz"))
end


# function components_on_plane(vf::PlanarVelocityField)
#     @unpack n1, n2, velocities, plane = vf
#     i, j = plane === :xy ? (1, 2) :
#            plane === :xz ? (1, 3) :
#            plane === :yz ? (2, 3) :
#            throw(ArgumentError("plane must be :xy, :xz, or :yz"))

#     U = reshape(getindex.(vs, i), n1, n2)
#     V = reshape(getindex.(vs, j), n1, n2)
#     U, V
# end

function viz(vf::PlanarVelocityField)
    points2 = [_proj2(vf.plane, x) for x in vf.points]
    vels2   = [_proj2(vf.plane, v) for v in vf.velocities]

    fig = Figure()
    (xl, yl) = _labels(vf.plane)
    ax = Axis(fig[1,1], aspect=DataAspect(), xlabel=xl, ylabel=yl)

    spacing = norm(points2[1] - points2[2]) 
    ar = arrows2d!(ax, points2, vels2;
        color=norm.(vels2),
        normalize=true,
        colormap=:ice,
        lengthscale=0.6f0*spacing
    )
    Colorbar(fig[2,1], ar, vertical=false, label=L"u")
    fig
end

function stream(vf::PlanarVelocityField; colorscale=identity)
    @unpack plane, a_range, b_range, c, velocities = vf
    ind = plane === :xy ? (1,2) : (plane === :xz ? (1,3) : (2,3))
    n1, n2 = length(a_range), length(b_range)
    U = reshape(getindex.(velocities, ind[1]), n1, n2)
    V = reshape(getindex.(velocities, ind[2]), n1, n2)

    itp_u = interpolate((a_range, b_range), U, Gridded(Linear()))
    itp_v = interpolate((a_range, b_range), V, Gridded(Linear()))

    # Continuous velocity function
    u_func(x, y) = Point2(itp_u(x, y), itp_v(x, y))

    vel_mag = sqrt.(U.^2 .+ V.^2)
    fig = Figure()
    (xl, yl) = _labels(plane)
    ax = Axis(fig[1,1], xlabel=xl, ylabel=yl, aspect=DataAspect())
    hm = heatmap!(ax, a_range, b_range, vel_mag, colorscale=colorscale)  # colorrange=(vmin, vmax))
    streamplot!(ax, u_func, 
        [minimum(a_range), maximum(a_range)], 
        [minimum(b_range), maximum(b_range)], 
        color = c -> :white, 
        arrow_size=10, 
        linewidth=1.5, 
        maxsteps=20000
    )
    Colorbar(fig[1, 2], hm, label="v")
    fig
end

viz(fv::FluidVelocity, 
    a_range::AbstractVector{T}, 
    b_range::AbstractVector{T}; 
    c::T=0.0, 
    plane=:xy) where {T <: Number} = viz(PlanarVelocityField(fv, a_range, b_range; c=c, plane=plane))

stream(fv::FluidVelocity, 
    a_range::AbstractVector{T}, 
    b_range::AbstractVector{T}; 
    c::T=0.0, 
    plane=:xy,
    colorscale=identity) where {T <: Number} = stream(
        PlanarVelocityField(fv, a_range, b_range; c=c, plane=plane), 
        colorscale=colorscale
    )

## Trajectories

function viz(traj::Trajectory; step=1)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    # b3 = cross.(traj.b1, traj.b2)
    lines!(ax, traj.x, color=:red)
    # spacing = norm(traj.x[1] - traj.x[2])
    # @info "" spacing
    # arrows3d!(ax, traj.x[1:step:end], traj.b1[1:step:end], normalize=true, lengthscale=step*spacing)
    # arrows3d!(ax, traj.x[1:step:end], traj.b2[1:step:end], normalize=true, lengthscale=step*spacing)
    fig
end

function viz_xy(traj::Trajectory)
    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())
    lines!(ax, traj.x, color=:red)
    fig
end



