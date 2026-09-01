using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FFTW
using StaticArrays


a = 4.5 # radius
α = π/6 # angle to connect the flagellum
f = 1.0 # frequency

N_body = 117   # body force points
Q_body = 993   # body quadrature points
N_flagellum = 53
Q_flagellum = 417

body = Part(EllipsoidBody(a, a, a), N_body, Q_body; eps=0.1)

# The flagellum model is a sum of standing waves
# θ(s) = A01*exp(iϕ01)*sin(π*s/2) + A11*exp(iϕ11)*sin(3π*s/2) + A21*exp(iϕ21)*sin(5π*s/2) + A31*exp(iϕ31)*sin(7π*s/2)
# θ(s,t) = real(exp(i*ω*t)*θ(s)) + exp(-i*ω*t)*conj(θ(s)))
#
# each standing wave mode sin(n+0.5)pi*s has a complex amplitude and phase so they can oscillate out of phase

f_model = PlanarStandingWaveFlagellum(
    70.,                              # L
    2π*f,                             # ω
    0.0,                              # C: static curvature
    SVector(0.18, 0.13, 0.46, 0.36), # A: mode amplitudes [A01, A11, A21, A31]
    SVector(0.26, -1.76, -0.07, 1.61) # ϕ: mode phases    [ϕ01, ϕ11, ϕ21, ϕ31]
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

flagellum = Part(
    f_model, N_flagellum, Q_flagellum,
    eps=0.1,
    location=[a*cos(α), a*sin(α), 0.]
)

pterosperma = MicroSwimmer([body, flagellum])


### Fluid field

x_points = range(-50.0, 70., 30)
y_points = range(-60.0, 60., 30)
traj_prob = SwimmingTrajectoryProblem(pterosperma, t_final=1.0, saveat=0.05)
solve_problem!(traj_prob)


# Disturbance field is in the lab frame, but shifted along with the swimming trajectory
vf = TimeAveragedDisturbanceField(traj_prob, x_points, y_points)
fig = stream(vf)


#### Visualise the trajectory

traj_prob = SwimmingTrajectoryProblem(pterosperma, t_final=10.0)
solve_problem!(traj_prob)
animate(traj_prob, step=1)
