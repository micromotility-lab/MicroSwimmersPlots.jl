function animate(prob::SwimmingTrajectoryProblem;
    wall=false, 
    limits=nothing, 
    step=5,
    filename=nothing,
    framerate=20,
    elevation=π/36, 
    azimuth=π/4
)
    traj = prob.traj
    ts = prob.traj.t
    swimmer = prob.swimming_problem.swimmer
    # b3 = reduce(hcat, cross.(traj.b1, traj.b2))
    
    t = Observable(1)
    T = @lift(traj.x[1:$t])

    if isnothing(limits)
        r_max = maximum(norm, swimmer.points.force_pts)
        xs = getindex.(traj.x, 1)
        ys = getindex.(traj.x, 2)
        zs = getindex.(traj.x, 3)
        limits = (
            (minimum(xs) - r_max, maximum(xs) + r_max),
            (minimum(ys) - r_max, maximum(ys) + r_max),
            (minimum(zs) - r_max, maximum(zs) + r_max)
        )
    end

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data, limits=limits, elevation=elevation, azimuth=azimuth)
    # ax = LScene(fig[1, 1])
    obs = viz!(ax, swimmer)
    lines!(ax, T, color=:red, linewidth=0.5)
    yield()
    display(fig)
    
    on(t) do i
        move_boundary!(swimmer, traj.x[i], traj.b1[i], traj.b2[i], ts[i])
        update_buffer!(obs[], swimmer)
        notify(obs)
    end

    # pts = @lift(get_points($s))
    # scatter!(ax, @lift(Point3f.(eachcol($s.config.quad_pts))), markersize=5)

    # lines!(ax, Point3f.(eachcol(X0)), color=:red, linewidth=0.5, linestyle=:dash)
    if wall
        v = [
            -1. -1. 0.;
            -1.  1. 0.;
            1.  1. 0.;
            1. -1. 0.
        ]
        f = [1 2 3; 3 4 1]
        mesh!(ax, v, f)
    end

    if isnothing(filename)
        for i in 1:step:length(ts)
            t[] = i
            # yield()
            sleep(1 // framerate)
        end
    else
        record(fig, filename, 1:length(ts); framerate=framerate) do i 
            t[] = i
        end 
    end
    # fig
end

"""Animate particle trajectories in a flow, NEEDS UPDATING"""
function animate(prob::ParticleTrajectoryProblem;
    wall=false, 
    limits=(-5., 5., -5., 5., -5., 5.), 
    step=5,
    filename=nothing,
    framerate=20,
    traj_length=600,
    elevation=π/6, 
    azimuth=π/4,
    num_particles=nothing
)
    traj = prob.trajectories
    ts = prob.t
    swimmer = prob.resistance_problem.swimmer

    t = Observable(1)
    num_particles = isnothing(num_particles) ? size(traj,1) ÷ 3 : num_particles

    particle_trajectories = [@lift(traj[3i-2:3i, max(1, $t-traj_length):$t]) for i in 1:num_particles]
    c = to_color(:purple)
    tailcol = [RGBAf(c.r, c.g, c.b, (i/traj_length)^2) for i in 1:traj_length]



    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data, limits=limits, elevation=elevation, azimuth=azimuth)
    # ax = LScene(fig[1, 1])
    obs = viz!(ax, swimmer)


    for i in 1:num_particles
        lines!(ax, particle_trajectories[i], linewidth=2.)
        scatter!(ax, @lift(Point3($(particle_trajectories[i])[:,end])), markersize=6)
    end

    if wall
        v = [
            -1. -1. 0.;
            -1.  1. 0.;
                1.  1. 0.;
                1. -1. 0.
        ]
        faces = [1 2 3; 3 4 1]
        mesh!(ax, v, faces)
    end
    display(fig)

    
    on(t) do i
        update_boundary!(swimmer, ts[i])
        update_buffer!(obs[], swimmer)
        notify(obs)
    end

    if isnothing(filename)
        for i in 1:step:length(ts)
            t[] = i
            sleep(0.05)
        end
    else
        record(fig, filename, 1:length(ts); framerate=framerate) do i 
            t[] = i
        end 
    end
end