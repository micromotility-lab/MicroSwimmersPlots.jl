using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FFTW


a = 4.5 # radius
α = π/6 # angle to connect the flagellum
f = 1.0 # frequency

N_body = 117   # body force points 
Q_body = 993   # body quadrature points
N_flagellum = 53 
Q_flagellum =417 

body = SphericalBody(a=4.5, N=N_body, Q=Q_body)

# The flagellum model is sum of standing waves
# θ(s) = A01*exp(iϕ01)*sin(π*s/2) + A11*exp(iϕ11)*sin(3π*s/2) + A21*exp(iϕ21)*sin(5π*s/2) + A31*exp(iϕ31)*sin(7π*s/2)
# θ(s,t) = real(exp(i*ω*t)*θ(s)) + exp(-i*ω*t)*conj(θ(s)))
# 
# each standing wave mode sin(n+0.5)pi*s has a complex amplitude and phase so they can oscillate out of phase

f_model = StandingWaveFlagellum(
    70.,    # L
    0.0,    # C
    0.18,   # A01
    0.26,   # ϕ01
    0.13,   # A11
    -1.76,  # ϕ01
    0.46,   # A21
    -0.07,  # ϕ21
    0.36,   # A31
    1.61,   # ϕ31
    2π*f    # ω
)

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

f = Flagellum(
    f_model, N_flagellum, Q_flagellum, 
    location=[a*cos(α), a*sin(α), 0.]
)

pterosperma = Flagellate(body, [f])


### Fluid field

x_points = range(-50.0, 70., 30)
y_points = range(-60.0, 60., 30)
traj_prob = SwimmingTrajectoryProblem(pterosperma, eps=0.1, t_final=1.0, saveat=0.05)
solve_problem!(traj_prob)


# Disturbance field is in the lab frame, but shifted along with the swimming trajectory
vf = TimeAveragedDisturbanceField(traj_prob, x_points, y_points)
fig = viz_xy(vf)


#### Visualise the trajectory

traj_prob = SwimmingTrajectoryProblem(pterosperma, t_final=10.0, eps=0.1)
solve_problem!(traj_prob)
animate(traj_prob, step=1)
