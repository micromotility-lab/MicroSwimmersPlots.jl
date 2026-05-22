using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using LinearAlgebra
# set_theme!(theme_dark())
set_theme!(theme_light())

#################################################################################
#### Example 1: An isolated flagellum ###########################################
#################################################################################

# Use a planar flagellum model for the tangent angle (Geyer 2016):
# θ(s,t) = Cs + (R₀ + R₁ sin(ks/L))*cos(ωt - ϕs/L)

model = PlanarFlagellum(
    1.,    #  L: Length of the flagellum (um)
    0.,    #  C: Static curvature 
    0.3,   #  R₀: Amplitude/envelope magnitude 
    0.15,  #  R₁: Spatial modulation of amplitude
    2π,    #  k: wavenumber of envelope modulation
    2π,    #  ϕ: wavenumber of travelling wave
    2π,    #  ω: angular frequency
    0.0    #  δ: overall phase (for multiple flagella)
)

# Next we construct a discretised Flagellum structure and apply our model

N = 23    # force points
Q = 127   # quadrature points Q > 4N

f = Flagellum(
    model,       # model: the model we defined above,
    N,             # N: number of force points
    Q              # Q: number of quadrature points
)

# Some alternative flagella beating envelopes changing R₀, R₁ and k

f = Flagellum(PlanarFlagellum(1., 0., 0.6, 0., 2π, 2π, 2π, 0.0), N, Q) 
f = Flagellum(PlanarFlagellum(1., 0., 1.3, 0., 2π, 2π, 2π, 0.0), N, Q)
f = Flagellum(PlanarFlagellum(1., 0., 0., 0.5, π/2, 2π, 2π, 0.0), N, Q)
f = Flagellum(PlanarFlagellum(1., 0., 0.6, 0.5, π/2, -2π, 2π, 0.0), N, Q)


# I've set things up so that you can run viz to visualise most objects

fig1 = viz(f)   # hold the left mouse button to move the camera around

# To solve the swimming problem (i.e. to calculate the rigid body velocity
# U, angular velocity Ω and force distribution) we set up a SwimmingProblem

prob = SwimmingProblem(f)

# This didn't solve the problem yet, we need to call solve_problem!

solve_problem!(prob)

# The ! is a julia convention which means that the argument you pass gets modified.
# We can get the results as follows

U = get_U(prob)
Ω = get_Ω(prob)
forces = get_forces(prob)

# Check that the total force and torque is zero 

F, T = total_force_and_torque(prob)  
 
# Plot the force distribution on the fluid due to the flagellum
fig2 = Figure()
ax = Axis(fig2[1,1], aspect=DataAspect())   # DataAspect draws equal axes
viz!(ax, f)                               # add the flagellum to the plot
ar = arrows2d!(ax, get_force_pts(prob), forces, color=norm.(forces), normalize=true, lengthscale=0.1)  # add quiver plot
arrows2d!(ax, Point3f(0.), U, color=:red)
Colorbar(fig2[2,1], ar, vertical=false, label="force")

# We can also get the velocity field at some predefined points in the x-y plane

x_points = range(-10.0, 11.0, 30)       # 30 equally spaced points between -10 and 11.0
y_points = range(-10.5, 10.5, 30) 

u = FluidVelocity(prob)
# stream plot
stream(u, x_points, y_points)
# arrow plot
vf = PlanarVelocityField(prob, x_points, y_points)
viz(vf)


# That was the velocity field at a particular time point. To see the average
# velocity field, construct a TimeAveragedVelocityField.

ave_vf = TimeAveragedPlanarVelocityField(
    prob, 
    x_points, 
    y_points; 
    period=1.0,       # Since we set ω=2π the period is 1.0 time unit
    num_t=30         # The number of time slices to include in the average
)

# Notice that the field is characteristic of a puller. What kinds of fields 
# do you get for the other flagella defined above?

fig3 = stream(ave_vf)

# Next let's see the swimming trajectory of the flagellum. To solve a time dependent
# swimming problem, construct a SwimmingTrajectoryProblem

tprob = SwimmingTrajectoryProblem(f, t_final=1.0, saveat=0.01) # solve for one period 100 pts per period
solve_problem!(tprob, periodic=true)

# The trajectory results are stored in tprob.traj
lines(tprob.traj.x)
viz(tprob)   # slide through the  trajectory



# For a periodic trajectory you can continue the trajectory by adding copies onto the end

traj = continue_periodic_trajectory(tprob.traj, 10) # continue for 10 periods

lines(traj.x)

# Finally watch an animation by calling animate on the SwimmingTrajectoryProblem
animate(traj, f)


##########################################################################################
#### Example 2. A sperm-like swimmer #####################################################
##########################################################################################

# Now we will create a cell body and attach a flagellum. We'll make a sphere of diameter 2um

a = 1.0 # semi-major axis 
b = 1.0 # semi-minor axis
c = 1.0 # semi-minor axis 

N_body = 213 # number of force points
Q_body = 917 # number of quadrature points

body = CellBody(EllipsoidBody(a, b, c), N, Q)

# Let's redefine the flagellum above, this time we will add a location where we want the flagellum to be attached

f = Flagellum(
    PlanarFlagellum(10., 0., 0.6, 0.5, π/2, 2π, 2π, 0.0),
    N,
    Q,
    location=[1.0, 0.0, 0.0],                         # connect to the edge on the x-axis 
    orientation=rotation_matrix([1.0, 0., 0.], 0.0)   # you can also change the orientation using rotation_matrix(axis, angle) 
)

# Put the body and flagellum together into a Flagellate structure

flg = Flagellate(
    body,
    [f]       # a list of flagella, with one element in this case
)

# average velocity field around the swimmer
sprob = SwimmingProblem(flg)
solve_problem!(sprob)
viz(sprob)

x_points = range(-10, 15,50)
y_points = range(-10, 10, 50)
ave_vf = TimeAveragedPlanarVelocityField(sprob, x_points, y_points)
stream(ave_vf)

tprob2 = SwimmingTrajectoryProblem(flg, t_final=1.0, saveat=0.05)
solve_problem!(tprob2, periodic=true)

# Again we continue the periodic trajectory to save unnecessary computation, and then animate
traj = continue_periodic_trajectory(tprob2.traj, 20)
animate(traj, flg)









