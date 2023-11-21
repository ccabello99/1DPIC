include("Grid.jl")
include("Fields.jl")
using LinearAlgebra
using Distributions

# Physical constants
kb = 1.38e-23   # Boltzmann constant
e = 1.6e-19     # Electron charge
me = 9.11e-31       # Electron mass
mp = 1.67e-27       # Proton mass
mi = 28 * mp        # Silicon ion mass
ϵ0 = 8.85e-12     # Vacuum permittivity

# Plasma parameters
#Te = 1e3*11604    # Initial temperatures (1 keV)
#Ti = 1e3*11604
nc = ϵ0 * me * ω0^2 / e^2   # Critical density (m^-3)
ni0 = 4.995e27              # Atomic density of Si
ne0 = 40*nc                 # Fully ionized Si electron plasma density

# Initialize dynamical variables
global ne = zeros(Nz)
global ni = zeros(Nz)
global ϵr = ones(Nz)
global ϵ = ϵ0.*ϵr
global ρ = zeros(Nz)
global Φ = zeros(Nz)
global Ez = zeros(Nz)
global Jy = zeros(Nz)
global Jz = zeros(Nz)

# Particle structure
mutable struct Particle
    z::Float64  # Position [z]
    v::Vector{Float64}  # Velocity [vz, vy]
    q::Float64  # Charge
    m::Float64  # Mass
    s::String   # Species
    w::Float64  # Weight
    γ::Float64  # Lorentz factor
end

# Initialize plasma density profile
function initGrad(ne, ni)

    L = 0.05 * λ
    front = 0.9 * Lz
    eplasma_gradient(x) = ne0 * exp((x - front) / L)
    iplasma_gradient(x) = ni0 * exp((x - front) / L)

    for i in 1:Nz
        if i*dz <= front
            ne[i] = eplasma_gradient(i * dz)+1e-36
            ni[i] = iplasma_gradient(i * dz)+1e-36
        else    
            ne[i] = ne0
            ni[i] = ni0 
        end
    end
    return ne, ni
end

# Initialize random velocities
function initVelocities(T, m)

    thermal_velocity = √(kb * T / m)
    dist = thermal_velocity * Normal(0, 1)
    vz = rand(dist)
    vy = 0
    v = [vz, vy]
    return v
end

function updateDensity(particles, ne)

    for particle in particles
        i = Int(floor(particle.z / dz)) 
        if i >= 1 && i < Nz
            ne[i] = particle.w / dz
        end
    end
    return ne
end

function assignWeight(particles)

    # Count number of macro particles in each grid cell
    local num_macroparticles = zeros(Int, Nz)

    for particle in particles
        i = Int(floor(particle.z / dz))  
        num_macroparticles[i] += 1    
    end

    # Assign weights based on macro particles per cell and species density
    for particle in particles
        i = Int(floor(particle.z / dz))

        if particle.s == "Electron"
            particle.w = ne[i] * dz / num_macroparticles[i]
        end

        if particle.s == "Ion"
            particle.w = ni[i] * dz / num_macroparticles[i]
        end

    end
    return particles
end

# Initialize particles
function createParticles(num_particles)

    particles = Particle[]
    electrons = Particle[]
    ions = Particle[]

    # Regions of (in)homogenous plasma density profile
    homog = Uniform(9/10, 1)
    inhomog = Uniform(6/10, 9/10)

    for i in 1:2num_particles
        # Initial velocity to start at rest
        v = [0, 0]
        γ = 1

            if i <= 0.75 * num_particles
                z = rand(homog) * Lz
                #v = initVelocities(Te, me)
                q = -e 
                w = 1.0
                m = me
                s = "Electron"
                push!(electrons, Particle(z, v, q, m, s, w, γ))
            
            elseif i <= 1.5 * num_particles
                z = rand(inhomog) * Lz
                #v = initVelocities(Te, me)
                q = -e
                w = 1.0
                m = me
                s = "Electron"
                push!(electrons, Particle(z, v, q, m, s, w, γ))

            elseif i <= 1.75 * num_particles
                z = rand(homog) * Lz
                #v = initVelocities(Ti, mi)
                q = 14*e
                w = 1.0
                m = mi
                s = "Ion"
                push!(ions, Particle(z, v, q, m, s, w, γ))

            else
                z = rand(inhomog) * Lz
                #v = initVelocities(Ti, mi)
                q = 14*e
                w = 1.0
                m = mi
                s = "Ion"
                push!(ions, Particle(z, v, q, m, s, w, γ))
            end 

        push!(particles, Particle(z, v, q, m, s, w, γ))
    end
    particles = assignWeight(particles)

    return particles, electrons, ions
end

# Linear interpolation to calculate charge/current density on grid
function InterpolateCharge(particles, ρ, Jz, Jy)

    for particle in particles
        i = Int(floor(particle.z / dz))
        if i >= 1 && i < Nz
            δz = particle.z / dz - i

            ρ[i] = abs(1 - δz) * (particle.q * particle.w / dz)
            ρ[i + 1] = abs(δz) * (particle.q * particle.w / dz)

            Jz[i] = ρ[i] * particle.v[1]
            Jz[i + 1] = ρ[i + 1] * particle.v[1]

            Jy[i] = ρ[i] * particle.v[2]
            Jy[i + 1] = ρ[i + 1] * particle.v[2]
        end
    end
    return ρ, Jz, Jy
end

# Poisson solver for 1D
function solvePoisson(ρ, ϵ)

    # Generate triagiagonal matrix for the linear system:
    dl = [-ones(Nz-2); 0]
    du = [0; -ones(Nz-2)]
    d = [1; 2*ones(Nz-2); 1]
    tri = Tridiagonal(dl, d, du)
    tri = tri ./ dz^2

    # Dirichlet BC
    ρ[1] = 0
    ρ[end] = maximum(ρ)

    # Solve
    rhs = ρ ./ ϵ
    Φ = tri \ rhs;

    return Φ
end

# Update static field
function updateStatic(ρ, Φ, ϵ, Ez)

    Φ = solvePoisson(ρ, ϵ)
    for i in 2:Nz-1
        Ez[i] -= (0.5 * dz * (Φ[i + 1] - Φ[i - 1]))
    end

    # Dirichlet BC
    Ez[1] = Ez[2]
    Ez[end] = Ez[end-1]

    return Φ, Ez
end

# Velocity update function
function velocityUpdate(particles, E_z, E_y, H_x)

    # Convert H to B for Lorentz force
    B_x =  μ0 .* H_x

    for particle in particles
        i = Int(floor(particle.z / dz))
        if i >= 1 && i < Nz && round(particle.w, digits=3) != 0
            δz = particle.z / dz - i

            # Linear interpolation of fields to grid
            Bx = (1 - δz) * B_x[i] + δz * B_x[i + 1]
            Ez = (1 - δz) * E_z[i] + δz * E_z[i + 1]
            Ey = (1 - δz) * E_y[i] + δz * E_y[i + 1]

            # Derivative of Ey^2 (assuming dy=dz) for ponderomotive force
            dEy2dy = 0.5 * (E_y[i + 1]^2 - E_y[i-1]^2) / dz


            # Velocity update with (Vay 2008) formalism

            # Initial values / create 3D vectors
            r = particle.q * dt / (particle.m)
            uz = particle.γ .* particle.v[1]
            uy = particle.γ .* particle.v[2]
            u = [0, uy, uz]
            v = [0, uy, uz] ./ particle.γ
            E = [0, Ey, Ez]
            B = [Bx, 0, 0]
            
            # Half of the velocity update 
            u_star = u + r * (E + 0.5 .* cross(v, B))

            # Compute auxilary values for final update
            tau = 0.5 * r .* B
            γprime = √(1 + (dot(u_star, u_star) ./ c0^2))
            σ = (γprime^2 - dot(tau, tau))
            w = dot(u_star, tau)
            particle.γ = √((σ + √(σ^2 + (dot(tau, tau) + w^2) ) ))
            t = tau ./ particle.γ
            s = 1 / (1 + dot(t, t))
        
            # Compute second-half of Vay update
            u = s * (u_star + dot(u_star, t) .* t + cross(u_star, t))

            # Include ponderomotive force (i.e. nonlinear term)
            u -= (dt * particle.q / (4 * particle.m * ω0^2)) .* [0, dEy2dy, 0]

            # Convert back to velocity
            particle.v[1] = u[3] ./ particle.γ
            particle.v[2] = (u[2] ./ particle.γ)

        end
    end
        
    return particles 

end

# Update particle positions and velocities
function updateParticles(particles, E_z, E_y, H_x) 

    # Perform half of position update
    for particle in particles
        i = Int(floor(particle.z / dz))
        if i >= 1 && i < Nz
            particle.z += 0.5 * particle.v[1] * dt
        end
    end

    # Perform velocity update
    particles = velocityUpdate(particles, E_z, E_y, H_x)

    # Perform second half of position update
    for particle in particles
        i = Int(floor(particle.z / dz))
        if i >= 1 && i < Nz
            particle.z += 0.5 * particle.v[1] * dt
        end
    end

    return particles
end

# Compute relative permittivity from plasma frequency
function updatePermittivity(ω, ne)
    ωp = zeros(Nz)
    
    for i in 1:Nz
        ωp[i] = sqrt(ne[i] * e^2 / (me * ϵ0))
    end

    ϵr = 1.0 .- ωp ./ ω

    return ϵr
end