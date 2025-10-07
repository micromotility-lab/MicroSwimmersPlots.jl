const I3f = Mat3f(I)

abstract type PlotBuffer end

function update_pts!(out_vec::Vector{Point3f}, in_vec::Vector{Point3f}, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    @inbounds for i in eachindex(in_vec)
        out_vec[i] = location + orientation*in_vec[i]
    end 
end

function update_pts!(out_vec::Vector{Point3f}, in_vec::Matrix{<:Number}, location::Point3f=Point3f(0.), orientation::Mat3f=I3f)
    @inbounds for i in axes(in_vec, 2)
        out_vec[i] = location + orientation*Point3f(in_vec[:,i])
    end 
end

function viz(fb::FluidBoundary)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    viz!(ax, fb)
    fig
end


# Cell bodies

struct CellBodyBuffer <: PlotBuffer
    mesh::GeometryBasics.Mesh
    ref_pts::Vector{Point3f}
end

function get_buffer(body::CellBody)
    mesh = gen_mesh(body)
    ref_pts = copy(coordinates(mesh))
    CellBodyBuffer(mesh, ref_pts)
end

function update_buffer!(buf::CellBodyBuffer, body::CellBody)
    update_pts!(buf.mesh.position, buf.ref_pts, Point3f(body.points.location), Mat3f(body.points.orientation))
end

function viz!(ax, buf::CellBodyBuffer)
    B = Observable(buf)
    mesh!(ax, @lift($B.mesh))
    B
end


function viz!(ax, body::CellBody)
    buf = get_buffer(body)
    update_buffer!(buf, body)
    viz!(ax, buf)
end


# Flagellum 

struct FlagellumBuffer <: PlotBuffer
    pts::Vector{Point3f}
end

get_buffer(f::Flagellum) = FlagellumBuffer(Vector{Point3f}(undef, f.points.N))

function update_buffer!(buf::FlagellumBuffer, f::Flagellum)
    update_pts!(buf.pts, f.points.force_pts, Point3f(f.points.location), Mat3f(f.points.orientation))
end

function update_buffer_observable!(B::Observable{FlagellumBuffer}, f::Flagellum)
    update_buffer!(B[], f)
    notify(B)
end

function viz!(ax, buf::FlagellumBuffer; linewidth=3, color=:forestgreen)
    B = Observable(buf)
    lines!(ax, @lift($B.pts), linewidth=linewidth, color=color)
    B
end

function viz!(ax, f::Flagellum)
    buf = get_buffer(f)
    update_pts!(buf.pts, f.points.force_pts, Point3f(f.points.location), Mat3f(f.points.orientation))
    viz!(ax, buf)
end



# Flagellate

struct FlagellateBuffer <: PlotBuffer
    body_buffer::CellBodyBuffer
    flagella_buffers::Vector{FlagellumBuffer}
end

function get_buffer(flg::Flagellate)
    body_buf = get_buffer(flg.body)
    flagella_bufs = [get_buffer(f) for f in flg.flagella]
    FlagellateBuffer(body_buf, flagella_bufs)
end

function update_buffer!(buf::FlagellateBuffer, flg::Flagellate)
    update_pts!(buf.body_buffer.mesh.position, buf.body_buffer.ref_pts, Point3f(flg.points.location), Mat3f(flg.points.orientation))
    for (i,f) in enumerate(flg.flagella)
        fb = buf.flagella_buffers[i]
        update_buffer!(fb, f)
        update_pts!(fb.pts, fb.pts, Point3f(flg.points.location), Mat3f(flg.points.orientation))
    end
end

function update_buffer_observable!(B::Observable{FlagellateBuffer}, flg::Flagellate)
    update_buffer!(B[], flg)
    notify(B)
end


function viz!(ax, buf::FlagellateBuffer, linewidth=3, color=:forestgreen)
    B = Observable(buf)
    mesh!(ax, @lift($B.body_buffer.mesh))
    for i in eachindex(buf.flagella_buffers)
        lines!(ax, @lift($B.flagella_buffers[i].pts), linewidth=linewidth, color=color)
    end
    B
end
  
function viz!(ax, flg::Flagellate)
    buf = get_buffer(flg)
    update_buffer!(buf, flg)
    viz!(ax, buf)
end


# Velocity fields

function viz!(ax, vf::VelocityField)
    spacing = norm(vf.points[1] - vf.points[2])
    arrows3d!(
        ax, 
        vf.points, 
        vf.velocities, 
        color=norm.(vf.velocities), 
        colormap=:ice,
        normalize=true, 
        lengthscale=spacing
    )
end


function viz(vf::VelocityField)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)  
    ar = viz!(ax, vf)
    Colorbar(fig[2,1], ar, vertical=false, label=L"u")
    fig
end

function viz_velocity_field_2d!(ax, points, velocities)
    spacing = norm(points[1] - points[2])
    arrows2d!(
        ax, 
        points, 
        velocities, 
        color=norm.(velocities), 
        normalize=true, 
        colormap=:ice,
        lengthscale=0.6*spacing
    )
end 

function viz_velocity_field_2d(points, velocities)
    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())
    ar = viz_velocity_field_2d!(ax, points, velocities)
    Colorbar(fig[2,1], ar, vertical=false, label=L"u")
    fig
end

function viz_xy(vf::VelocityField)
    points = [Point2f(x[1], x[2]) for x in vf.points]
    velocities = [Point2f(v[1], v[2]) for v in vf.velocities]
    viz_velocity_field_2d(points, velocities)
end

function viz_xz(vf::VelocityField)
    points = [Point2f(x[1], x[3]) for x in vf.points]
    velocities = [Point2f(v[1], v[3]) for v in vf.velocities]
    viz_velocity_field_2d(points, velocities)
end

function viz_yz(vf::VelocityField)
    points = [Point2f(x[2], x[3]) for x in vf.points]
    velocities = [Point2f(v[2], v[3]) for v in vf.velocities]
    viz_velocity_field_2d(points, velocities)
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


