# 1D PIC
using Random
using Plots
using Distributions
using Colors
using LinearAlgebra
using FFTW
using SparseArrays

include("Grid.jl")
include("Particles.jl")
include("Fields.jl")

# Initialize plasma
num_particles = 1000
global particles, electrons, protons = createParticles(num_particles)

# Calculate Charge Distribution, Static Fields, and Interaction Parameters
global ρ, Jz, Jy = InterpolateCharge(particles, ρ, Jz, Jy)

global Φ, Ez = updateStatic(ρ, Φ, ϵ, Ez)

global ϵr = updatePermittivity(ρ, ω0, ne)

#display(plot(z, Φ.*8, label="Electric Potential"))
#display(scatter!(z, ρ, label="Charge Density"))
#end

#T = range(1, nsteps, nsteps)
#display(plot(T, ge.(T), xlims=[0, nsteps]))
#end


#eyr = ones(Complex{Float64}, nfreq)
#eyt = ones(Complex{Float64}, nfreq)
#src = ones(Complex{Float64}, nfreq)

for t = 1:nsteps

    # Calculate permittivity of the grid
    global ϵr = updatePermittivity(ρ, ω0, ne)

    # Solve Maxwell's eqs
    global Ey, Hx, E1, E2, H1, H2 = FDTD(Ey, Hx, E1, E2, H1, H2, Jy, t)

    particle_z = []
    #for particle in particles
    #    push!(particle_z, particle.z)
    #end
    #print(particle_z)
    # Particle Push
    global particles = updateParticles(particles, Ez, Ey, Hx) 

    # Calculate Charge/Current Density and Static Fields
    global ρ, Jz, Jy = InterpolateCharge(particles, ρ, Jz, Jy)

    global Φ, Ez = updateStatic(ρ, Φ, ϵ, Ez)


    # Plotting
    if t % 5 == 0
        #p = plot(z, Ey / E0, legend = true, title = "1D PIC Simulation", lw = 2, xlabel = "Z (μm)", label = "Ey")
        #plot!(z, Hx, lw = 2, label = "Hx")
        #p = plot(z, Ez, label="Ez (static)")
        #ρ
        p = scatter(z, ρ, label="Charge Density")
        #scatter!(z, Jz, label="Current Density (z)")
        #scatter!(z, Jy, label="Current Density (y)")
        #plot!(z, Φ)
        # Draw dielectric
        #plot!([1050, 1050, 2000, 2000], [-1.5, 1.5, 1.5, -1.5], lw = 3, color = :black, legend = false)
        display(p)
        sleep(0.02)
    end
    #display()

    # Draw dielectric
    #display()


end

