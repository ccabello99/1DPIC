include("Grid.jl")
include("Fields.jl")
using LinearAlgebra

# Physical constants
kb = 1.38e-23   # Boltzmann constant
#e = 1.6e-19     # Electron charge
#me = 9.11e-31       # Electron mass
#mp = 1.67e-27       # Proton mass
# Atomic units
me = 1
mp = 1800
e = 1

# Plasma parameters
Te_i = 200    # Initial temperatures
Tp_i = 200
#ne0 = 198 * 1.74e21   # Species densities (m^-3)
ne0 = 6.6e29
global ne = zeros(Nz)
ni0 = ne0
global ni = zeros(Nz)
global ϵr = ones(Nz)
ϵ0 = 8.85e-12     # Vacuum permittivity
global ϵ = ϵ0.*ϵr
#cs = √(kb * Te_i / mp)  # Ion sound velocity
#global L(t) = cs * t           # Plasma density expansion 

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

# Initialize plasma density gradient
L = 0.5 * λ
front = 0.9 * Lz
eplasma_gradient(x) = ne0 * exp((x - front) / L)
iplasma_gradient(x) = ni0 * exp((x - front) / L)

for i in 1:Nz
    if i*dz <= front
        ne[i] = eplasma_gradient(i * dz)
        ni[i] = iplasma_gradient(i * dz)
    else
        ne[i] = ne0
        ni[i] = ni0 
    end
end

# Initialize vectors
global ρ = zeros(Nz)
global Φ = zeros(Nz)
global Ez = zeros(Nz)
global Jy = zeros(Nz)
global Jz = zeros(Nz)

function assignWeight(particles)

    # First, count number of macro particles in each grid cell
    local num_macroparticles = zeros(Int, Nz)

    for particle in particles    
        # Determine the cell indices
        i = Int(floor(particle.z / dz))  

        # Increase the count of macro-particles in the cell
        num_macroparticles[i] += 1    
    end

    # Now, assign weights based on number of macro particles per cell
    for particle in particles

        # Get macro_particle indices
        i = Int(floor(particle.z / dz))

        if particle.s == "Electron"
            particle.w = ne[i] * dz / num_macroparticles[i]
        end

        if particle.s == "Proton"
            particle.w = ni[i] * dz / num_macroparticles[i]
        end

end
    return particles
end

# Initialize particles
function createParticles(num_particles)
    particles = Particle[]
    electrons = Particle[]
    protons = Particle[]
    homog = Uniform(9/10, 1)
    inhomog = Uniform(5/10, 9/10)

    for i in 1:2num_particles
        # Initial position
        #z = rand(homog) * Lz
        #z = (0.25+rand()*0.01)* Lz     # For testing
        # Initial velocity
        v = [0, 0]
        γ = 1

            if i <= 0.75 * num_particles
                z = rand(homog) * Lz
                #z = (0.05+rand()*0.25)* Lz
                q = -e 
                w = 1.0
                m = me
                s = "Electron"
                push!(electrons, Particle(z, v, q, m, s, w, γ))
            
            elseif i <= 1.5 * num_particles
                z = rand(inhomog) * Lz
                #z = (0.6+rand()*0.15)* Lz
                q = -e
                w = 1.0
                m = me
                s = "Electron"
                push!(electrons, Particle(z, v, q, m, s, w, γ))

            elseif i <= 1.75 * num_particles
                z = rand(homog) * Lz
                #z = (0.6+rand()*0.15)* Lz
                q = e
                w = 1.0
                m = mp
                s = "Proton"
                push!(protons, Particle(z, v, q, m, s, w, γ))

            else
                z = rand(inhomog) * Lz
                q = e
                w = 1.0
                m = mp
                s = "Proton"
                push!(protons, Particle(z, v, q, m, s, w, γ))
            end 

        push!(particles, Particle(z, v, q, m, s, w, γ))
    end
    #particles = assignWeight(particles)

    return particles, electrons, protons
end

#print(ne)
#print(dz)
#display(scatter(z, ne))
#print(ne[2])

function InterpolateCharge(particles, ρ, Jz, Jy)

    for particle in particles
        i = Int(floor(particle.z / dz))
        #i = Int(clamp(round(particle.z / dz), 1, Nz))
        if i >= 1 && i < Nz && round(particle.w, digits=3) != 0
            δz = particle.z / dz - i

            ρ[i] = abs(1 - δz) * (particle.q * particle.w / dz)
            ρ[i + 1] = abs(δz) * (particle.q * particle.w / dz)

            Jz[i] += ρ[i] * particle.v[1]
            Jz[i + 1] += ρ[i + 1] * particle.v[1]

            Jy[i] += ρ[i] * particle.v[2]
            Jy[i + 1] += ρ[i + 1] * particle.v[2]
        end
    end
    return ρ, Jz, Jy
end

# Define the Poisson solver function
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

    # Solve:
    rhs = ρ ./ ϵ
    Φ = tri \ rhs;

    return Φ
end

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

function borisUpdate(particle, E_z, E_y, H_x)
    # Convert H to B for Lorentz force
    B_x =  μ0 .* H_x

    for particle in particles
        i = Int(floor(particle.z / dz))
        #i = Int(clamp(round(particle.z / dz), 1, Nz))
        if i >= 1 && i < Nz && round(particle.w, digits=3) != 0
            δz = particle.z / dz - i

            # Interpolate fields to grid
            Bx = (1 - δz) * B_x[i] + δz * B_x[i + 1]
            Ez = (1 - δz) * E_z[i] + δz * E_z[i + 1]
            Ey = (1 - δz) * E_y[i] + δz * E_y[i + 1]

            # Update coefficient
            qmdt = 0.5 * particle.q * dt / (particle.m)

            # # Define 3D vector for relativistic momentum
            uz = particle.v[1] * particle.γ
            uy = particle.v[2] .* particle.γ
            ux = 0
            u3 = [0, uy, uz]
    
            # First half of velocity update
            u_minus = u3 + qmdt * [0, Ey, Ez]
    

            # Auxilary values for rotation
            particle.γ = √(1 + (dot(u_minus,u_minus) / c0^2))
            t = qmdt .* [Bx, 0, 0] ./ particle.γ
            s = 2 * t / (1 + dot(t, t))

            # Velocity rotation plus second half of update
            u_prime = u_minus + cross(u_minus, t)
            u_plus = u_minus + cross(u_prime, s)

            # Final update of each particle's velocity
            particle.v[1] = u_plus[3] + qmdt * Ez
            particle.v[2] = u_plus[2] + qmdt * Ey
        end
    end
    return particles
end

function velocityUpdate(particles, E_z, E_y, H_x)

    # Convert H to B for Lorentz force
    B_x =  μ0 .* H_x

    for particle in particles
        i = Int(floor(particle.z / dz))
        if i >= 1 && i < Nz && round(particle.w, digits=3) != 0
            δz = particle.z / dz - i

            # Interpolate fields to grid
            Bx = (1 - δz) * B_x[i] + δz * B_x[i + 1]
            Ez = (1 - δz) * E_z[i] + δz * E_z[i + 1]
            Ey = (1 - δz) * E_y[i] + δz * E_y[i + 1]


            # Velocity update with (Vay 2008) formalism

            # Initial values
            #γ = 1 / √(1 - (norm(particle.v) / c0^2))
            u = particle.γ .* particle.v
            r = particle.q * dt / (2*particle.m)

            # Half of the velocity update 
            u_half = u + r * ([Ez, Ey] + 0.5 .* [-(u[2]/particle.γ)*Bx, (u[1]/particle.γ)*Bx])
            
            # Compute other half (u_star zero since Bx . [Ez, Ey] = 0)
            u_prime = u_half +  r .* [Ez, Ey]
            u_star = 0

            # Compute auxilary values for final update
            tau = r * Bx
            γprime = √(1 + (norm(u_prime) / c0^2) )
            σ = (γprime^2 - tau^2)
            particle.γ = √((σ + √(σ^2 + 4*(tau^2 + u_star^2) ) ) / 2)
            t = tau ./ particle.γ
            s = 1 / (1 + t^2)
        
            # Final update
            u = s * (u_prime + [-u_prime[2] * t, u_prime[1] * t])
            particle.v = u ./ particle.γ 
        end
    end
        
    return particles 

end

function updateParticles(particles, E_z, E_y, H_x) 

    # Perform half position update
    for particle in particles
        i = Int(floor(particle.z / dz))
        #i = Int(clamp(round(particle.z / dz), 1, Nz))
        if i >= 1 && i < Nz
            particle.z += 0.5 * particle.v[1] * dt
        end
    end

    # Perform velocity update
    #particles = velocityUpdate(particles, E_z, E_y, H_x)
    particles = borisUpdate(particles, E_z, E_y, H_x)

    # Perform second half of position update
    for particle in particles
        i = Int(floor(particle.z / dz))
        #i = Int(clamp(round(particle.z / dz), 1, Nz))
        if i >= 1 && i < Nz
            particle.z += 0.5 * particle.v[1] * dt
        end
    end

    return particles
end

function updatePermittivity(ρ, ω, ne)
    ωp = zeros(Nz)
    
    for i in 1:Nz
        ωp[i] = sqrt(ne[i] * e^2 / (me * ϵ0))
    end

    ϵr = 1.0 .- ωp ./ ω

    return ϵr
end