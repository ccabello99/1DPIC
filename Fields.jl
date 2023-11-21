include("Grid.jl")

# Physical Constants
ϵ0 = 8.85e-12     # Vacuum permittivity
μ0 = 4*π*1e-7     # Vacuum permeability
η0 = √(μ0/ϵ0)     # Vacuum impedence
c0 = 1 / √(ϵ0*μ0)
me = 9.11e-31
e = 1.6e-19

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
global H3 = H2
global E1 = 0
global E2 = E1
global E3 = E2

# Update coefficients
S = 1/2         # Courant Number
loss = zeros(Nz)
global σ = zeros(Nz)
global ceh = zeros(Nz)
global che = zeros(Nz)
global cee = zeros(Nz)

# Update FDTD coeffiecients (includes loss term for plasma)
function updateCoefficients(ϵr, μr)
    for i in 1:Nz
        σ[i] = (ne[i] * e^2) / me
        loss[i] = 0.5 * σ[i] * dt / ϵ[i]
        cee[i] = (1.0 - loss[i]) / (1.0 + loss[i]) 
        ceh[i] = c0 * dt ./ ϵr[i] ./ dz ./ (1.0 + loss[i])
        che[i] = c0 * dt ./ μr[i] ./ dz
    end
    return ceh, che, cee, σ
end

# Position of the source
z_init = Int(Nz / abs(Nz) + 5)

# Initialize source parameters
ω0 = 2π * 374.74057 * 1e12 # 800 nm
k = (ω0 / c0)
λ = 2π / k
Nλ = λ / dz
f = (2π * c0 * dt) / (Nλ * dz)
#τ = 36  # Pulse duration (~FWHM)
τ = 512
#τ = 640
t0 = 6 * τ


# Charge density parameters
nc = ϵ0 * me * ω0^2 / e^2   # Critical density (m^-3)
ni0 = 4.995e27              # Atomic density of Si
ne0 = 40*nc              # Fully ionized Si electron plasma density


## Gaussian Envelope Parameters
# Normalized relativistic field amplitude
a0 = 0.001      
# Electric field amplitude
E0 = (c0 * me * ω0 * a0) / e              
# H field time/spatial grid offset
n_src = sqrt(ϵr[z_init] * μr[z_init])
delta_t = (n_src * dz) / (2 * c0) + (dt / 2)    
# Gaussian pulse centered at ω0 with correction to ΔCEP=0
ge(t) = E0 .* exp(-((t - t0) ./ (τ)).^2) * cos((f * t)+π*5/4)
gh(t) = -(E0 / η0) .* exp(-((t - t0 + delta_t) ./ (τ)).^2) * cos((f * t)+π*5/4)


# Ricker Wavelet
arg = S / Nλ
Md = 20
fre(t) = E0 .* (1 - 2*(π^2)*(arg * (t - Md))^2) * exp(-π^2 * (arg * (t - Md))^2)
frh(t) = A .* (1 - 2*(π^2)*(arg * (t - Md + delta_t))^2) * exp(-π^2 * (arg * (t - Md + delta_t))^2)


#FDTD Loop for one time step
function FDTD(Ey, Hx, E1, E2, E3, H1, H2, H3, Jy, t)

    # H field loop
    for k = 1:Nz-1
        Hx[k] += che[k] * (Ey[k+1] - Ey[k])
    end

    # Absorbing Boundary Condition at end of grid
    Hx[Nz] += che[Nz] * (E3 - Ey[Nz])

    # H field source
    Hx[z_init-1] -= che[z_init-1] * gh(t)

    # Record H at boundary
    H3, H2, H1 = H2, H1, Hx[1]
    
    # E field loop
    for k = 2:Nz
        Ey[k] = cee[k] * Ey[k] + ceh[k] * ((Hx[k] - Hx[k-1]))
    end

    #Source current correction
    Ey = Ey - ((dt ./ (ne .* e^2)) .* Jy)

    # Absorbing Boundary Condition at start of grid
    Ey[1] += ceh[1] * ((Hx[1] - H3))

    # Record E at boundary
    E3, E2, E1 = E2, E1, Ey[Nz]

    # E Field source
    Ey[z_init] -= ceh[z_init] * ge(t)

    return Ey, Hx, E1, E2, E3, H1, H2, H3
end