include("Grid.jl")
#include("Particles.jl")

# Physical Constants
ϵ0 = 8.85e-12     # Vacuum permittivity
μ0 = 4*π*1e-7     # Vacuum permeability
c0 = 1 / √(ϵ0*μ0)
#me = 9.11e-31
#e = 1.6e-19
# Atomic units
me = 1
e = 1

# Dynamical variables
global Ey = zeros(Nz)
global Hx = zeros(Nz)
global ϵr = ones(Nz)
global μr = ones(Nz)
global ϵ = ϵ0.*ϵr

# Uncomment to set ϵr and μr for specific regions
#ϵr[1050:2000] .= 4.2
#μr[1050:2000] .= 0.9

# To record boundaries
global H1 = 0
global H2 = H1
global E1 = 0
global E2 = E1

# Update coefficients
η0 = √(μ0/ϵ0)
S = 1/2
function updateCoefficients(ϵr, μr)
    ce = c0 * dt ./ ϵr ./ dz
    ch = c0 * dt ./ μr ./ dz
    return ce, ch
end


# Position of the source
z_init = Int(Nz / abs(Nz) + 4)

# Initialize source

ω0 = 2π * 374.74057 * 1e12 # 800 nm
k = (ω0 / c0)
λ = 2π / k
Nλ = λ / dz
f = (2π * c0 * dt) / (Nλ * dz)
τ = 32
#τ = 640
t0 = 3 * τ

a0 = 1e-4      # Normalized field amplitude
E0 = (c0 * me * ω0 * a0) / e              # Electric field amplitude
A = -E0 .* sqrt(ϵr[z_init] / μr[z_init])      # Magnetic field amplitude
ge(t) = E0 .* exp(-((t - t0) ./ (τ)).^2) * cos(f * t)
n_src = sqrt(ϵr[z_init] * μr[z_init])
delta_t = (n_src * dz) / (2 * c0) + (dt / 2)
gh(t) = A .* exp(-((t - t0 + delta_t) ./ (τ)).^2) * cos(f * t)

#FDTD Loop for one time step
function FDTD(Ey, Hx, E1, E2, H1, H2, Jy, t)
    
    # Set update coefficients
    ce, ch = updateCoefficients(ϵr, μr)

    # Record H at boundary
    H2, H1 = H1, Hx[1]

    # H field loop
    for k = 1:Nz-1
        Hx[k] += ch[k] * (Ey[k+1] - Ey[k])
    end

    Hx[Nz] += ch[Nz] * (E2 - Ey[Nz])

    # Record E at boundary
    E2, E1 = E1, Ey[Nz]

    Ey[1] += ce[1] * (Hx[1] - H2)

    # E field loop
    for k = 2:Nz
        Ey[k] += ce[k] * ((Hx[k] - Hx[k-1]) - Jy[k])
    end

    # E Field source
    Ey[z_init] -= ce[z_init] * gh(t)

    # H field source
    Hx[z_init-1] -= ch[z_init-1] * ge(t)

    return Ey, Hx, E1, E2, H1, H2
end