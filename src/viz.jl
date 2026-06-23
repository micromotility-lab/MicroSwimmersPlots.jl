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

function update_pts!(out_vec::Vector{Point3f}, in_vec::AbstractVector{<:AbstractVector}, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    @inbounds for i in eachindex(in_vec)
        out_vec[i] = location + orientation * Point3f(in_vec[i])
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

# for a cell body, this does not apply its local location and orientation
function viz!(ax, fb::FluidBoundary, location=nothing, orientation=nothing; kwargs...)
    buf = get_buffer(fb)
    if !isnothing(location) && !isnothing(orientation)
        update_buffer!(buf, fb, Point3f(location), Mat3f(orientation))
    else
        update_buffer!(buf, fb)
    end
    viz!(ax, buf; kwargs...)
end


# Cell bodies

mutable struct CellBodyBuffer <: PlotBuffer
    mesh::GeometryBasics.Mesh
    ref_pts::Vector{Point3f}
end

function get_buffer(part::Part{<:CellBodyModel}; slice=false)
    mesh = slice ? gen_mesh_sliced(part.model) : gen_mesh(part.model)
    ref_pts = copy(coordinates(mesh))
    CellBodyBuffer(mesh, ref_pts)
end

function update_buffer!(buf::CellBodyBuffer, part::Part{<:CellBodyModel}, location=Point3f(0.), orientation=I3f)
    update_pts!(buf.mesh.position, buf.ref_pts, Point3f(part.frame.location), Mat3f(part.frame.orientation))
    update_pts!(buf.mesh.position, buf.mesh.position, location, orientation)
end

function update_mesh!(buf::CellBodyBuffer, part::Part{<:CellBodyModel}, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    buf.mesh = gen_mesh(part.model)
    buf.ref_pts = copy(coordinates(buf.mesh))
    update_pts!(buf.mesh.position, buf.ref_pts, Point3f(part.frame.location), Mat3f(part.frame.orientation))
    update_pts!(buf.mesh.position, buf.mesh.position, location, orientation)
end

viz!(ax, B::Observable{CellBodyBuffer}; kwargs...) = mesh!(ax, @lift($B.mesh); kwargs...)

# Flagellum 

struct BareFlagellumBuffer <: FlagellumBuffer
    pts::Vector{Point3f}
end

get_buffer(f::Part{<:FlagellumModel}) = BareFlagellumBuffer(Vector{Point3f}(undef, nq(f.disc)))

function update_buffer!(buf::BareFlagellumBuffer, f::Part{<:FlagellumModel}, location=Point3f(0.), orientation=I3f)
    update_pts!(buf.pts, f.disc.quad_pts, Point3f(f.frame.location), Mat3f(f.frame.orientation))
    update_pts!(buf.pts, buf.pts, Point3f(location), Mat3f(orientation))
end

function viz!(ax, B::Observable{BareFlagellumBuffer}; linewidth=3, color=:forestgreen)
    lines!(ax, @lift($B.pts), linewidth=linewidth, color=color)
end

# MicroSwimmer

struct MicroSwimmerBuffer <: PlotBuffer
    body_buffers::Vector{CellBodyBuffer}
    flagellum_buffers::Vector{BareFlagellumBuffer}
end

function get_buffer(ms::MicroSwimmer)
    body_bufs = [get_buffer(p) for p in ms.parts if p.model isa CellBodyModel]
    flag_bufs = [get_buffer(p) for p in ms.parts if p.model isa FlagellumModel]
    MicroSwimmerBuffer(body_bufs, flag_bufs)
end

function update_buffer!(buf::MicroSwimmerBuffer, ms::MicroSwimmer,
    location=Point3f(ms.frame.location),
    orientation=Mat3f(ms.frame.orientation))
    body_parts = filter(p -> p.model isa CellBodyModel, ms.parts)
    flag_parts = filter(p -> p.model isa FlagellumModel, ms.parts)
    for (i, p) in enumerate(body_parts)
        update_buffer!(buf.body_buffers[i], p, location, orientation)
    end
    for (i, p) in enumerate(flag_parts)
        update_buffer!(buf.flagellum_buffers[i], p, location, orientation)
    end
end

function viz!(ax, B::Observable{MicroSwimmerBuffer};
    bodycolor=Makie.wong_colors()[1],
    linewidth=3,
    color=:forestgreen,
    rasterize_body=true
)
    for i in eachindex(B[].body_buffers)
        viz!(ax, @lift($B.body_buffers[i]); color=bodycolor, rasterize=rasterize_body)
    end
    for i in eachindex(B[].flagellum_buffers)
        viz!(ax, @lift($B.flagellum_buffers[i]); linewidth=linewidth, color=color)
    end
end

# problems

function viz!(ax, prob::SwimmingProblem; rasterize=1, step=1, kwargs...)
    B = viz!(ax, prob.microswimmer; rasterize_body=rasterize, kwargs...)
    N_body = sum(nf(p.disc) for p in prob.microswimmer.parts if p.model isa CellBodyModel; init=0)
    forces = get_forces(prob)[N_body+1:step:end]
    force_pts = get_force_pts(prob)[N_body+1:step:end]
    Fmag = norm.(forces)
    Fmax = maximum(Fmag)
    flag_idx = findfirst(p -> p.model isa FlagellumModel, prob.microswimmer.parts)
    L = isnothing(flag_idx) ? 1.0 : prob.microswimmer.parts[flag_idx].model.L

    ar = arrows2d!(
        ax,
        force_pts,
        forces,
        color = (:darkred,0.8),                  # or :darkred
        normalize = true,
        lengthscale = 0.25 * L .* (Fmag ./ Fmax),
        rasterize=rasterize
    )

    # cs = [(Makie.wong_colors()[2], norm(f)/M) for f in forces]
    # ar = arrows2d!(ax, get_force_pts(prob), forces, color=cs, normalize=true, lengthscale=0.1*L)  # add quiver plot
    B
end

viz(prob::SwimmingProblem; kwargs...) = begin
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    viz!(ax, prob; kwargs...)
    fig
end

# function viz(prob::SwimmingProblem)
#     ax = Axis3(fig[1,1], aspect=:data)
#     B = viz!(ax, prob.microswimmer)
#     forces = get_forces(prob)
#     M = maximum(norm.(forces))
#     cs = [(Makie.wong_colors()[2], norm(f)/M) for f in forces]
#     L = prob.microswimmer.flagella[1].model.L
    
#     ar = arrows2d!(ax, get_force_pts(prob), forces, color=cs, normalize=true, lengthscale=0.1*L)  # add quiver plot
#     # Colorbar(fig[2,1], ar, vertical=false, label="force")
#     display(fig)
#     B
# end

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
    plane === :xy && return (L"x \, (\mu\mathrm{m})", L"y \, (\mu\mathrm{m})")
    plane === :xz && return (L"x \, (\mu\mathrm{m})", L"z \, (\mu\mathrm{m})")
    plane === :yz && return (L"y \, (\mu\mathrm{m})", L"z \, (\mu\mathrm{m})")
    throw(ArgumentError("plane must be :xy, :xz, or :yz"))
end

function viz(vf::PlanarVelocityField; fig = Figure())
    points2 = [_proj2(vf.plane, x) for x in vf.points]
    vels2   = [_proj2(vf.plane, v) for v in vf.velocities]

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

# function stream!(parent, vf::PlanarVelocityField; 
#     show_colorbar=true, 
#     stream_kwargs=(;arrow_size=10, linewidth=1.5, rasterize=1),
#     hm_kwargs=(;rasterize=1),

# )
#     g = parent isa GridLayout ? parent : GridLayout(parent)
#     @unpack plane, a_range, b_range, c, velocities = vf
#     ind = plane === :xy ? (1,2) : (plane === :xz ? (1,3) : (2,3))
#     n1, n2 = length(a_range), length(b_range)
#     U = reshape(getindex.(velocities, ind[1]), n1, n2)
#     V = reshape(getindex.(velocities, ind[2]), n1, n2)

#     itp_u = interpolate((a_range, b_range), U, Gridded(Linear()))
#     itp_v = interpolate((a_range, b_range), V, Gridded(Linear()))

#     # Continuous velocity function
#     u_func(x, y) = Point2(itp_u(x, y), itp_v(x, y))

#     # vel_mag = sqrt.(U.^2 .+ V.^2)
#     # vel_mag = [norm(u_func(x,y)) for x in range(a_range[1], a_range[end], 300), y in range(b_range[1], b_range[end], 300)]
   
    # vel_mag_raw = sqrt.(U.^2 .+ V.^2)  # compute at original grid points
    # itp_mag = interpolate((a_range, b_range), vel_mag_raw, Gridded(Linear()))
    # a_itp = range(a_range[1], a_range[end], 300)
    # b_itp = range(b_range[1], b_range[end], 300)
    # vel_mag = [itp_mag(x, y) for x in a_itp, y in b_itp]
        
#     (xl, yl) = _labels(plane)
#     ax = Axis(g[1,1], xlabel=xl, ylabel=yl, aspect=DataAspect())
#     hm = heatmap!(ax, a_itp, b_itp, vel_mag; hm_kwargs...)  # colorrange=(vmin, vmax))
#     streamplot!(ax, u_func, 
#         [minimum(a_range), maximum(a_range)], 
#         [minimum(b_range), maximum(b_range)]; 
#         color = c -> :white, 
#         maxsteps=20000,
#         stream_kwargs...
#     )
#     if show_colorbar
#         Colorbar(g[1, 2], hm, width=5)
#         a = (b_range[end] - b_range[1]) / (a_range[end] - a_range[1])
#         rowsize!(g, 1, Aspect(1, a))
#         # Label(g[1,2,Top()], L"v \, (\mu\mathrm{m/beat})")
#         Label(g[1,2,Top()], L"v \, (\mu\mathrm{m/sec})")
#         colgap!(g, 5)
#     end
#     # colsize!(fig.layout, 1, Aspect(1, 1.0))
#     # rowsize!(fig.layout, 1, Aspect(1, 1.0))
#     ax
# end

function stream!(parent, vf::PlanarVelocityField;
    show_colorbar=true,
    in_domain=(x, z) -> true,
    stream_kwargs=(; arrow_size=10, linewidth=1.5, rasterize=1),
    hm_kwargs=(; rasterize=1),
)
    g = parent isa GridLayout ? parent : GridLayout(parent)
    @unpack plane, a_range, b_range, c, velocities = vf
    ind = plane === :xy ? (1,2) : (plane === :xz ? (1,3) : (2,3))
    n1, n2 = length(a_range), length(b_range)
    U = reshape(getindex.(velocities, ind[1]), n1, n2)
    V = reshape(getindex.(velocities, ind[2]), n1, n2)

    # Mask outside domain

    vel_mag_raw = sqrt.(U.^2 .+ V.^2)  # compute at original grid points
    itp_mag = interpolate((a_range, b_range), vel_mag_raw, Gridded(Linear()))
    a_itp = range(a_range[1], a_range[end], 300)
    b_itp = range(b_range[1], b_range[end], 300)
    vel_mag = [itp_mag(x, y) for x in a_itp, y in b_itp]


    # vel_mag = sqrt.(U.^2 .+ V.^2)
    for (i, x) in enumerate(a_itp), (j, z) in enumerate(b_itp)
        if !in_domain(x, z)
            vel_mag[i, j] = NaN
        end
    end

    itp_u = interpolate((a_range, b_range), U, Gridded(Linear()))
    itp_v = interpolate((a_range, b_range), V, Gridded(Linear()))

    u_func(x, z) = begin
        !in_domain(x, z) && return Point2(NaN, NaN)
        u, v = itp_u(x, z), itp_v(x, z)
        sqrt(u^2 + v^2) < 1e-6 && return Point2(NaN, NaN)
        Point2(u, v)
    end

    # u_func(x, z) = in_domain(x, z) ?
    #     Point2(itp_u(x, z), itp_v(x, z)) :
    #     Point2(0.0, 0.0)

    (xl, yl) = _labels(plane)
    ax = Axis(g[1,1], xlabel=xl, ylabel=yl, aspect=DataAspect())
    hm = heatmap!(ax, a_itp, b_itp, vel_mag;
        hm_kwargs...
    )
    streamplot!(ax, u_func,
        [minimum(a_range), maximum(a_range)],
        [minimum(b_range), maximum(b_range)];
        color=c -> :white,
        maxsteps=20000,
        stream_kwargs...
    )
    if show_colorbar
        Colorbar(g[1,2], hm, width=5)
        a = (b_range[end] - b_range[1]) / (a_range[end] - a_range[1])
        rowsize!(g, 1, Aspect(1, a))
        Label(g[1,2,Top()], L"v \, (\mu\mathrm{m/sec})")
        colgap!(g, 5)
    end
    ax
end

stream(vf::PlanarVelocityField; kwargs...) = begin
    fig = Figure()
    stream!(fig[1,1], vf; kwargs...)
    fig
end

viz(fv::FluidVelocity, 
    a_range::AbstractVector{T}, 
    b_range::AbstractVector{T}; 
    c::T=0.0, 
    plane=:xy) where {T <: Number} = viz(PlanarVelocityField(fv, a_range, b_range; c=c, plane=plane))

stream!(parent, fv::FluidVelocity, 
    a_range::AbstractVector{T}, 
    b_range::AbstractVector{T}; 
    c::T=0.0, 
    plane=:xy,
    kwargs...) where {T <: Number} = stream!(parent,
        PlanarVelocityField(fv, a_range, b_range; c=c, plane=plane); 
        kwargs...
)

stream(fv::FluidVelocity, 
    a_range::AbstractVector{T}, 
    b_range::AbstractVector{T}; 
    c::T=0.0, 
    plane=:xy,
    kwargs...) where {T <: Number} = stream(
        PlanarVelocityField(fv, a_range, b_range; c=c, plane=plane); 
        kwargs...
)

# function vorticity!(parent, vf::PlanarVelocityField;
#     show_colorbar=true,
#     method=:interpolant,  # :interpolant or :finitediff
#     colorrange=nothing,
#     show_streamlines=true,
#     stream_kwargs=(; arrow_size=10, linewidth=1.5, rasterize=1),
#     hm_kwargs=(; rasterize=1),
# )
#     g = parent isa GridLayout ? parent : GridLayout(parent)
#     @unpack plane, a_range, b_range, c, velocities = vf
#     ind = plane === :xy ? (1,2) : (plane === :xz ? (1,3) : (2,3))
#     n1, n2 = length(a_range), length(b_range)
#     U = reshape(getindex.(velocities, ind[1]), n1, n2)
#     V = reshape(getindex.(velocities, ind[2]), n1, n2)

#     itp_u = interpolate((a_range, b_range), U, Gridded(Linear()))
#     itp_v = interpolate((a_range, b_range), V, Gridded(Linear()))

#     ω = if method === :interpolant
#         _vorticity_interpolant(itp_u, itp_v, a_range, b_range)
#     else
#         _vorticity_finitediff(U, V, a_range, b_range)
#     end

#     ω_max = isnothing(colorrange) ? maximum(abs, ω) : colorrange
#     (xl, yl) = _labels(plane)
#     ax = Axis(g[1,1], xlabel=xl, ylabel=yl, aspect=DataAspect())
#     hm = heatmap!(ax, a_range, b_range, ω;
#         colormap=:RdBu,
#         colorrange=(-ω_max, ω_max),
#         hm_kwargs...
#     )

#     if show_streamlines
#         u_func(x, y) = Point2(itp_u(x, y), itp_v(x, y))
#         streamplot!(ax, u_func,
#             [minimum(a_range), maximum(a_range)],
#             [minimum(b_range), maximum(b_range)];
#             color=c -> :black,
#             maxsteps=20000,
#             stream_kwargs...
#         )
#     end

#     if show_colorbar
#         Colorbar(g[1,2], hm, width=5)
#         a = (b_range[end] - b_range[1]) / (a_range[end] - a_range[1])
#         rowsize!(g, 1, Aspect(1, a))
#         Label(g[1,2,Top()], L"\omega \, (\mathrm{s}^{-1})")
#         colgap!(g, 5)
#     end
#     ax
# end



function vorticity!(parent, vf::PlanarVelocityField;
    show_colorbar=true,
    in_domain=(x, z) -> true,
    method=:interpolant,
    colorrange=nothing,
    show_streamlines=true,
    stream_kwargs=(; arrow_size=10, linewidth=1.5, rasterize=1),
    hm_kwargs=(; rasterize=1),
)
    g = parent isa GridLayout ? parent : GridLayout(parent)
    @unpack plane, a_range, b_range, c, velocities = vf
    ind = plane === :xy ? (1,2) : (plane === :xz ? (1,3) : (2,3))
    n1, n2 = length(a_range), length(b_range)
    U = reshape(getindex.(velocities, ind[1]), n1, n2)
    V = reshape(getindex.(velocities, ind[2]), n1, n2)

    itp_u = interpolate((a_range, b_range), U, Gridded(Linear()))
    itp_v = interpolate((a_range, b_range), V, Gridded(Linear()))

    a_itp = range(a_range[1], a_range[end], 300)
    b_itp = range(b_range[1], b_range[end], 300)

    ω = if method === :interpolant
        _vorticity_interpolant(itp_u, itp_v, a_itp, b_itp)
    else
        _vorticity_finitediff(U, V, a_itp, b_itp)
    end


    # Mask outside domain
    for (i, x) in enumerate(a_itp), (j, z) in enumerate(b_itp)
        if !in_domain(x, z)
            ω[i, j] = NaN
        end
    end

    ω_max = isnothing(colorrange) ? maximum(abs, filter(!isnan, ω)) : colorrange
    (xl, yl) = _labels(plane)
    ax = Axis(g[1,1], xlabel=xl, ylabel=yl, aspect=DataAspect())
    hm = heatmap!(ax, a_itp, b_itp, ω;
        colormap=:RdBu,
        colorrange=(-ω_max, ω_max),
        nan_color=:transparent,
        hm_kwargs...
    )

    if show_streamlines
        u_func(x, z) = in_domain(x, z) ?
            Point2(itp_u(x, z), itp_v(x, z)) :
            Point2(0.0, 0.0)
        streamplot!(ax, u_func,
            [minimum(a_range), maximum(a_range)],
            [minimum(b_range), maximum(b_range)];
            color=c -> :black,
            maxsteps=20000,
            stream_kwargs...
        )
    end

    if show_colorbar
        Colorbar(g[1,2], hm, width=5)
        a = (b_range[end] - b_range[1]) / (a_range[end] - a_range[1])
        rowsize!(g, 1, Aspect(1, a))
        Label(g[1,2,Top()], L"\omega \, (\mathrm{s}^{-1})")
        colgap!(g, 5)
    end
    ax
end

vorticity(vf::PlanarVelocityField; kwargs...) = begin
    fig = Figure()
    vorticity!(fig[1,1], vf; kwargs...)
    fig
end

function _vorticity_interpolant(itp_u, itp_v, a_range, b_range)
    [
        Interpolations.gradient(itp_v, x, z)[1] - Interpolations.gradient(itp_u, x, z)[2]
        for x in a_range, z in b_range
    ]
end

function _vorticity_finitediff(U, V, a_range, b_range)
    dx = diff(a_range)
    dz = diff(b_range)
    dVdx = similar(V)
    dUdz = similar(U)

    # dV/dx
    for j in axes(V, 2)
        dVdx[1, j]   = (V[2, j]   - V[1, j])   / dx[1]
        dVdx[end, j] = (V[end, j] - V[end-1, j]) / dx[end]
        for i in 2:size(V, 1)-1
            dVdx[i, j] = (V[i+1, j] - V[i-1, j]) / (dx[i-1] + dx[i])
        end
    end

    # dU/dz
    for i in axes(U, 1)
        dUdz[i, 1]   = (U[i, 2]   - U[i, 1])   / dz[1]
        dUdz[i, end] = (U[i, end] - U[i, end-1]) / dz[end]
        for j in 2:size(U, 2)-1
            dUdz[i, j] = (U[i, j+1] - U[i, j-1]) / (dz[j-1] + dz[j])
        end
    end

    dVdx .- dUdz
end
    

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



