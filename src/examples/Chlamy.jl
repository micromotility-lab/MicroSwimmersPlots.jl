using MicroSwimmers
using MicroSwimmersPlots
using GLMakie


a=5.    # major axis length in um
b=4.    # minor axis length
c=4.    # minor axis length
L=10.   # Length of flagella
C=-2.5  # Intrinsic curvature of flagellum
α=π/16   # attachment angle of flagella in x-y plane relative to x-axis
N_body=411   # number of force points in body discretisation
Q_body=3111   # number of quadrature points in body discretisation
N_f=31       # number of force points in flagellum discretisation
Q_f=151       # number of quadrature points in flagellum discretisation

body = EllipsoidBody(a, b, c, N_body, Q_body)

# Use a planar flagellum model for the tangent angle:
# θ(s,t) = Cs + (R₀ + R₁ sin(ks/L))*cos(ωt - ϕs/L)
model = PlanarFlagellum(
    L,    #  L: length 
    C,    #  C: curvature
    0.3,  #  R₀: amplitude
    0.15, #  R₁: spatial modulation of amplitude 
    2π,   #  k: wavenumber of spatial modulation
    2π,   #  ω: angular frequency
    2π,   #  ϕ: wavenumber of travelling wave
    0.0   #  δ: overall phase (if you want asynchrony)
)

model = PlanarFlagellum(11., C, .8, 0.15, 1.0π, 2π, 2π, 0.0)


# 3D flagella with stokeslets arranged on a tube

f1_tube = TubeFlagellum(
    model,    # which model are you using for the flagellum? (PlanarFlagellum in this case)
    N_f,      # number of force points
    5,
    Q_f,      # number of quadrature points
    17,
    location=[a*cos(α), b*sin(α), 0],              # attachment point of flagellum base
    orientation=rotation_matrix([0.,0.,1.], π/4)*rotation_matrix([1.,0.,0.], 1.0π)  # orientation relative to reference config of the model (change the angle for tilt)
)

f2_tube = TubeFlagellum(
    model, N_f, 5, Q_f, 17,
    location=[a*cos(α), -b*sin(α), 0],
    orientation=rotation_matrix([0., 0., 1.], -π/4)*rotation_matrix([1.,0.,0.], 0.0),
)

# Centerline stokeslets
f1_line = Flagellum(
    model,    # which model are you using for the flagellum? (PlanarFlagellum in this case)
    N_f,      # number of force points
    Q_f,      # number of quadrature points
    location=[a*cos(α), b*sin(α), 0],              # attachment point of flagellum base
    orientation=rotation_matrix([0.,0.,1.], π/4)*rotation_matrix([1.,0.,0.], 4π/5)  # orientation relative to reference config of the model (change the angle for tilt)
)

f2_line = Flagellum(
    model, N_f, Q_f,
    location=[a*cos(α), -b*sin(α), 0],
    orientation=rotation_matrix([0., 0., 1.], -π/4)*rotation_matrix([1.,0.,0.], 0.0),
)


# Force points are centerline stokeslets, quadrature points on a tube
f1_line_tube = LineTubeFlagellum(
    model,    # which model are you using for the flagellum? (PlanarFlagellum in this case)
    N_f,      # number of force points
    Q_f,      # number of quadrature points
    17,
    location=[a*cos(α), b*sin(α), 0],              # attachment point of flagellum base
    orientation=rotation_matrix([1.,0.,0.], 1.0π)  # orientation relative to reference config of the model (change the angle for tilt)
)

f2_line_tube = LineTubeFlagellum(
    model, N_f, Q_f, 17,
    location=[a*cos(α), -b*sin(α), 0],
    orientation=rotation_matrix([1.,0.,0.], 0.0),
)



# Create a list containing the two flagella
flagella = [f1_line, f2_line]

# Finally return a flagellate with a body and list of flagella
chlamy = Flagellate(
    body,
    flagella
)


viz(chlamy)

prob = SwimmingProblem(chlamy)
x_points = range(-20.0, 21.0, 30)       # 30 equally spaced points between -10 and 11.0
y_points = range(-19.5, 19.5, 30) 

vf = TimeAveragedVelocityField(prob, x_points, y_points, 0.0)
fig_vf = viz_xy(vf)

# We'll solve for the swimming trajectory
prob = SwimmingTrajectoryProblem(chlamy, t_final=1.0, saveat=0.125, eps=1e-3)
solve_problem!(prob, periodic=true)

# Since the solution is periodic we can just replicate the trajectory 20 times
continue_periodic_trajectory!(prob.traj, 1000)
prob.traj.x[end]
viz(prob.traj)

animate(prob, step=13, azimuth=π/4, elevation=π/4) #, limits=(-5, 35, -10, 10, -8, 8))
animate(prob, step=13, azimuth=π/4, elevation=π/4, limits=(-5, 35, -10, 10, -8, 8), filename="chlamy_swimming.mp4", framerate=30) 


