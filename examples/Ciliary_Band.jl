using MicroSwimmers
using MicroSwimmersPlots
using GLMakie

# In this example we will use the same simple model for a tintinnid and platynereis.
# First we define a function which generates the microswimmer

function BandedSphericalMicroSwimmer(;
    a=4.,                    # radius of sphere
    body_N=109,              # force points
    body_Q=447,              # quadrature points
    num_flagella=32,
    total_phase_diff=12π,    # metachronal wavelength (12pi = 6 wavelengths)
    tilt_angle=0.,           # tilt of ciliary beat plane
    h=0.                     # the height of the band above/below the equator
)
    body = Part(EllipsoidBody(a, a, a), body_N, body_Q)

    # some maths
    R_tilt = rotation_matrix([1., 0., 0.], tilt_angle)
    R = rotation_matrix([0., 0., 1.0], π/2)
    r = sqrt(a^2 - h^2)
    θ = atan(h, r)

    # generate flagella positioned around the sphere
    flagella = [Part(
        PlanarFlagellum(1., -2.5, 0.7, 0.15, 2π, 2π, 2π, n*total_phase_diff/num_flagella),
        15,
        117,
        location=[h, r*cos(2π*n/num_flagella), r*sin(2π*n/num_flagella)],
        orientation = rotation_matrix([1., 0., 0.], 2pi*n/num_flagella) * rotation_matrix([0., 0., 1.], -θ) * R * R_tilt
    ) for n in 1:num_flagella]

    MicroSwimmer([body; flagella])
end

#############################################################################
#### Helical Swimming: Platynereis ##########################################
#############################################################################

# The default parameters above will give you a platynereis with 12 flagella
# and a metachronal wave with 6 wavelengths

platy = BandedSphericalMicroSwimmer()
viz(platy)

platy_prob = SwimmingTrajectoryProblem(platy, t_final=6.0, saveat=0.1)
solve_problem!(platy_prob, periodic=true)
continue_periodic_trajectory!(platy_prob.traj, 40)
animate(platy_prob)

# try changing the tilt angle to see the effect on the trajectory, it will run
# faster with fewer flagella

#############################################################################
#### Filter Feeder: Tintinnid ###############################################
#############################################################################

# Changing a few of the parameters we can produce something more like a mouth

tintinnid = BandedSphericalMicroSwimmer(
    a=2.,
    num_flagella=24,
    total_phase_diff=12π,
    tilt_angle=π/8,
    h=-1.7
)
viz(tintinnid)
tintinnid_prob = ParticleTrajectoryProblem(
    tintinnid,
    ys=range(-3,3,6),
    zs=range(-3,3,6),
    x=-4.5,
    t_final=20.0
)

solve_problem!(tintinnid_prob)
animate(tintinnid_prob)
