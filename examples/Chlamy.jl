using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FFTW


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

body = Part(EllipsoidBody(a, b, c), N_body, Q_body)

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


## You can use this function to get the spectrum for the standing wave flagellum from data
## shear_angle is a matrix containing one period of oscillation of size S x T
## with space points in the first dimension and time in the second
## leave off the initial point at s=0 (where shear_angle=0 anyway)
## The modes 01 - 31 above correspond to amplitude and phase of spectrum[1:4,2]
## try plotting spectrum with a heatmap to see if shear_angle is well captured by these modes

function spectrum_from_data(shear_angle)
    n, m = size(shear_angle)
    return FFTW.fft(FFTW.r2r(shear_angle, FFTW.RODFT01, 1), 2) / (n*m)
end


f1 = Part(
    model,
    N_f,
    Q_f,
    location=[a*cos(α), b*sin(α), 0],
    orientation=rotation_matrix([0.,0.,1.], π/4)*rotation_matrix([1.,0.,0.], 4π/5)
)

f2 = Part(
    model,
    N_f,
    Q_f,
    location=[a*cos(α), -b*sin(α), 0],
    orientation=rotation_matrix([0., 0., 1.], -π/4)*rotation_matrix([1.,0.,0.], 0.0),
)


chlamy = MicroSwimmer([body, f1, f2])


viz(chlamy)

prob = SwimmingProblem(chlamy)
x_points = range(-20.0, 21.0, 30)       # 30 equally spaced points between -10 and 11.0
y_points = range(-19.5, 19.5, 30)

ave_vf = TimeAveragedPlanarVelocityField(prob, x_points, y_points; period=1.0, num_t=30)
stream(ave_vf)

# We'll solve for the swimming trajectory
prob = SwimmingTrajectoryProblem(chlamy, t_final=1.0, saveat=0.125, eps=1e-3)
solve_problem!(prob, periodic=true)

# Since the solution is periodic we can just replicate the trajectory 50 times
continue_periodic_trajectory!(prob.traj, 50)
prob.traj.x[end]
viz(prob.traj)

animate(prob, step=13, azimuth=π/4, elevation=π/4)
# animate(prob, step=13, azimuth=π/4, elevation=π/4, limits=(-5, 35, -10, 10, -8, 8), filename="chlamy_swimming.mp4", framerate=30)
