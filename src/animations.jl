function animate(microswimmer::MicroSwimmer, T=5.0, num_t=20*5+1)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    B = viz!(ax, microswimmer)
    display(fig)
    for t in range(0, T, num_t)
        update_boundary!(microswimmer, t)
        update_buffer_observable!(B, microswimmer)
        sleep(0.01)
    end
end

function animate(
    traj::Trajectory,
    microswimmer::MicroSwimmer;
    wall=false, 
    limits=nothing, 
    step=5,
    filename=nothing,
    framerate=30,
    compression=20,
    elevation=π/36, 
    azimuth=π/4,
    fig_size=(1920, 1080)
)
    ts = traj.t
    
    t = Observable(1)
    T = @lift(traj.x[1:$t])

    if isnothing(limits)
        r_max = maximum(norm, microswimmer.points.force_pts)
        xs = getindex.(traj.x, 1)
        ys = getindex.(traj.x, 2)
        zs = getindex.(traj.x, 3)
        limits = (
            (minimum(xs) - r_max, maximum(xs) + r_max),
            (minimum(ys) - r_max, maximum(ys) + r_max),
            (minimum(zs) - r_max, maximum(zs) + r_max)
        )
    end

    fig = Figure(size=fig_size)
    # ax = LScene(fig[1, 1])
    ax = Axis3(fig[1, 1], aspect=:data, limits=limits, elevation=elevation, azimuth=azimuth)
    obs = viz!(ax, microswimmer)
    lines!(ax, T, color=:red, linewidth=0.5)
    yield()
    display(fig)
    
    on(t) do i
        move_boundary!(microswimmer, traj.x[i], traj.b1[i], traj.b2[i], ts[i])
        update_buffer_observable!(obs, microswimmer)
    end

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
        record(fig, filename, 1:length(ts); framerate=framerate, compression=compression) do i 
            t[] = i
        end 
    end
    # fig
end

animate(prob::SwimmingTrajectoryProblem; kwargs...) = animate(
        prob.traj, prob.swimming_problem.microswimmer; 
        kwargs...   
)


"""Animate particle trajectories in a flow, NEEDS UPDATING"""
function animate(
    ts::Vector,
    traj::Matrix,
    microswimmer::MicroSwimmer;
    wall=false, 
    limits=(-5., 5., -5., 5., -5., 5.), 
    step=5,
    traj_length=600,
    filename=nothing,
    framerate=20,
    elevation=π/6, 
    azimuth=π/4,
    num_particles=nothing
)
    t = Observable(1)
    num_particles = isnothing(num_particles) ? size(traj,1) ÷ 3 : num_particles

    particle_trajectories = [@lift(traj[3i-2:3i, max(1, $t-traj_length):$t]) for i in 1:num_particles]

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data, limits=limits, elevation=elevation, azimuth=azimuth)
    obs = viz!(ax, microswimmer)


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
        update_boundary!(microswimmer, ts[i])
        update_buffer_observable!(obs, microswimmer)
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

animate(prob::ParticleTrajectoryProblem; kwargs...) = animate(
        prob.t, prob.trajectories, prob.resistance_problem.microswimmer; 
        kwargs...   
)