# 1D PIC
using Random
using Plots
using Distributions
using Colors
using LinearAlgebra
using FFTW
using DelimitedFiles

include("Grid.jl")
include("Particles.jl")
include("Fields.jl")

# Initialize plasma
num_particles = 2000
global ne, ni = initGrad(ne, ni)
global particles, electrons, ions = createParticles(num_particles)

# Calculate Initial Charge Distribution and Static Fields
global ρ, Jz, Jy = InterpolateCharge(particles, ρ, Jz, Jy)
global Φ, Ez = updateStatic(ρ, Φ, ϵ, Ez)

# To visualize before time loops 
#display(plot(z, Φ, label="Electric Potential"))
#display(scatter!(z, ρ, label="Charge Density"))
#end

# Vector to save field for later processing
Ey_save = zeros(nsteps)
node = Int(round(0.5 * Nz))

for t = 1:nsteps

    # Update electron density and update permittivity
    global ne = updateDensity(particles, ne)
    global ϵr = updatePermittivity(ω0, ne)

    # Update FDTD coefficients
    global ceh, che, cee, σ = updateCoefficients(ϵr, μr)

    # Solve Maxwell's eqs
    global Ey, Hx, E1, E2, E3, H1, H2, H3 = FDTD(Ey, Hx, E1, E2, E3, H1, H2, H3, Jy, t)

    # Particle Push
    global particles = updateParticles(particles, Ez, Ey, Hx) 

    # Interpolate Charge/Current Density to grid
    global ρ, Jz, Jy = InterpolateCharge(particles, ρ, Jz, Jy)

    # Calculate Static Fields
    global Φ, Ez = updateStatic(ρ, Φ, ϵ, Ez)

    # Save Electric field component of interest for later processing 
    Ey_save[t] = Ey[node]

    # Plotting
    if t % 10 == 0
        local p = plot(z, Ey, label = "Ey", title = "1D PIC Simulation", lw = 2, ylims=[-1.5*E0, 1.5*E0], xlabel = "Z (μm)", ylabel="Field Strength (V/m)", legend = true)
        plot!(z, Hx, lw = 2, label = "Hx")
        #p = plot(z, Ez, label="Ez (static)")
        #p = scatter(particle_z, Ey)
        #scatter!(z, ρ ./ norm(ρ), label="Charge Density")
        #plot!(z, Jz, label="Current Density (z)")
        plot!(z, (dt ./ (ne .* e^2)) .* Jy, label = "Current Density (y)", color=:green)
        #plot!(z, Φ)
        display(p)
        #savefig(p, "plot_$t.png")
        sleep(0.025)
    end

end


# Save incident and reflected fields from Ey_save
#Ey_inc = Ey_save[1:1000]
#Ey_ref = Ey_save[1001:end]

#writedlm( "Ey_inc.csv",  Ey_inc, ',')
#writedlm( "Ey_ref.csv",  Ey_ref, ',')
