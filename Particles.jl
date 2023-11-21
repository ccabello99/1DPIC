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
Z = 14              # Silicon atomic number
ϵ0 = 8.85e-12     # Vacuum permittivity

# Plasma parameters
num_particles = 45      # Number of particles per cell
Te = 1e3*11604    # Initial temperatures (1 keV)
Ti = 1e3*11604
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

    L = 0.5 * λ
    front = 0.9 * Lz
    back = 1 * Lz
    eplasma_gradient(x) = ne0 * exp((x - front) / L)
    iplasma_gradient(x) = ni0 * exp((x - front) / L)

    for i in 1:Nz
        if i*dz <= front
            ne[i] = eplasma_gradient(i * dz)+1e-40
            ni[i] = iplasma_gradient(i * dz)+1e-40
        end
        if i*dz > front && i*dz <= back
            ne[i] = ne0
            ni[i] = ni0 
        end
    end
    return ne, ni
end

# Initialize random velocities
function initVelocities(T, m, q)

    thermal_velocity = √(abs(q/e) * kb * T / m)
    dist = thermal_velocity * Normal(0, 1)
    vz = rand(dist)
    vy = 0
    v = [vz, vy]
    return v
end

# Count number of macro particles in each grid cell
function countParticles(particles)

    local num_macroparticles = zeros(Int, Nz)

    for particle in particles
        i = Int(round(particle.z / dz))  
        num_macroparticles[i] += 1    
    end

    return num_macroparticles
end

# Update electron density
function updateDensity(particles, ne)

    num_macroparticles = countParticles(particles)

    for particle in particles
        i = Int(round(particle.z / dz)) 
        if i >= 1 && i <= Nz
            if particle.s == "Electron"
                ne[i] = particle.w * num_macroparticles[i] / dz
            end
        end
    end
    return ne
end

# Assign weights based on macro particles per cell and species density
function assignWeight(particles)

    num_macroparticles = countParticles(particles)

    for particle in particles
        i = Int(round(particle.z / dz))

        if particle.s == "Electron"
            particle.w = ne[i] * dz / num_macroparticles[i]
        end

        if particle.s == "Ion"
            particle.w = ni[i] * dz / num_macroparticles[i]
        end

    end
    return particles
end

# Create particles at a certain position
function createParticles(num_particles, pos)

    electrons = Particle[]
    ions = Particle[]

    # Regions of (in)homogenous plasma density profile
    #homog = Uniform(9/10, 1)
    #inhomog = Uniform(7/10, 9/10)

    for i in 1:num_particles
        # Initial velocity to start at rest
        #v = [0, 0]
        #γ = 1        

            if i <= 36
                z =  pos * Lz
                q = -e
                v = initVelocities(Te, me, q)
                w = 1.0
                m = me
                s = "Electron"
                γ = √(1 - dot(v, v) / c0^2)
                push!(electrons, Particle(z, v, q, m, s, w, γ))
            else
                z = pos * Lz
                q = Z*e
                v = initVelocities(Ti, mi, q)
                w = 1.0
                m = mi
                s = "Ion"
                γ = √(1 - dot(v, v) / c0^2)
                push!(ions, Particle(z, v, q, m, s, w, γ))
            end 
    end
    
    return electrons, ions
end

# Initialize all particles 
function initParticles()

    particles = Particle[]

    # Initial position corresponding to wanted grid index
    # Currently set for  particles (36 e & 9 ions) per cell
    pos = 0.7998046875
    last = 1
    #pos = 0.6500244140625
    #last = 0.85
    step = 1/8192
    num_iter = Int(round((last - pos) / step)) + 1
    
    for i in 1:num_iter
        electrons, ions = createParticles(num_particles, pos)
        # Append electrons and ions to full particles list
        for electron in electrons
            push!(particles, electron)
        end
        for ion in ions
            push!(particles, ion)
        end
        pos += step
    end

    particles = assignWeight(particles)

    return particles
end

# Linear interpolation to calculate charge/current density on grid
function InterpolateCharge(particles, ρ, Jz, Jy)

    for particle in particles
        i = Int(round(particle.z / dz))
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
        i = Int(round(particle.z / dz))
        if i >= 1 && i < Nz && round(particle.w, digits=3) != 0
            δz = particle.z / dz - i

            # Linear interpolation of fields to grid
            Bx = (1 - δz) * B_x[i] + δz * B_x[i + 1]
            Ez = (1 - δz) * E_z[i] + δz * E_z[i + 1]
            Ey = (1 - δz) * E_y[i] + δz * E_y[i + 1]

            # Central finite difference for ∇(Ey^2) (assuming dy=dz)
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
        i = Int(round(particle.z / dz))
        if i >= 1 && i < Nz
            particle.z += 0.5 * particle.v[1] * dt
        end
    end

    # Perform velocity update
    particles = velocityUpdate(particles, E_z, E_y, H_x)

    # Perform second half of position update
    for particle in particles
        i = Int(round(particle.z / dz))
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