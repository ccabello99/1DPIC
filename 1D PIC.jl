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
global ne, ni = initGrad(ne, ni)
global particles = initParticles()

# Calculate Initial Charge Distribution and Static Fields
global ρ, Jz, Jy = InterpolateCharge(particles, ρ, Jz, Jy)
global Φ, Ez = updateStatic(ρ, Φ, Ez)


## To visualize before time loops 
#num_macro = countParticles(particles)
#display(scatter(z, num_macro, label="Number of Macro Particles"))
#tt = range(1, nsteps, nsteps)
#display(plot(tt .* dt .* 1e15, ge.(tt) ./ E0, xlims=[0,200]))
#display(plot!(tt .* dt .* 1e15, gh.(tt) ./ (E0 / η0), xlims=[0,200]))
#display(plot(z, Φ, label="Electric Potential"))
#display(plot(z, ρ, label="Charge Density"))
#display(plot(z, ne, label="Electron Density"))
#display(plot!(z, ni, label="Ion Density"))
#print(length(particles))
#end

# Vector to save field for later processing
Ey_save = zeros(nsteps)
Jy_save = zeros(nsteps)
node = Int(round(0.5 * Nz))

for t = 1:nsteps

    # Update electron density and update permittivity
    global ne = updateDensity(particles, ne)
    global ϵr = updatePermittivity(ω0, ne)

    # Update FDTD coefficients
    global ceh, che, cee, σ = updateCoefficients(ne, ϵr, μr)

    # Solve Maxwell's eqs
    global Ey, Hx, E1, E2, E3, H1, H2, H3 = FDTD(Ey, Hx, E1, E2, E3, H1, H2, H3, Jy, t)

    # Particle Push
    global particles = updateParticles(particles, Ez, Ey, Hx) 

    # Interpolate Charge/Current Density to grid
    global ρ, Jz, Jy = InterpolateCharge(particles, ρ, Jz, Jy)

    # Calculate Static Fields
    global Φ, Ez = updateStatic(ρ, Φ, Ez)

    # Save Electric field component of interest for later processing 
    Ey_save[t] = Ey[node]

    # Plotting
    if t % 20 == 0
        local p = plot(z, 2 .* Ey ./ E0, label = "Ey", title = "1D PIC Simulation", lw = 2, xlabel = "Z (μm)", legend = true)
        ylims!(-1.5,1.5)
        #ylims!(-1.5*E0, 1.5*E0)
        plot!(z, 2 .* Hx ./ E0, lw = 2, label = "Hx")
        #p = plot(z, Ez, label="Ez")
        #p = scatter(particle_z, Ey)
        #scatter!(z, ρ ./ norm(ρ), label="Charge Density")
        plot!(z, ne ./ ne0, label="Electron Density", lw = 2, color=:red)
        #plot!(z, ni ./ maximum(ni), label="Ion Density")
        #local p = plot(z, (1 .- loss) ./ (1 .+ loss), label="Relative Permittivity")
        plot!(z, Jy ./ maximum(abs.(Jy)), label = "Jy", color=:green, lw = 2)
        hline!([0], color=:black, lw = 2.5, label="")
        #plot!(z, Φ)
        display(p)
        savefig(p, "plot_$t.png")
        #sleep(0.025) 
    end

end


# Save incident and reflected fields from Ey_save
half = Int(nsteps/2)
Ey_inc = Ey_save[1:half]
Ey_ref = Ey_save[half+1:end]

writedlm( "Ey_inc.csv",  Ey_inc, ',')
writedlm( "Ey_ref.csv",  Ey_ref, ',')
