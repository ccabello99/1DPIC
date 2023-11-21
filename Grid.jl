# Size of the grid
Nz = 8192
#Nz = 512       # For testing
Lz = 4e-5

# Grid parameters
c0 = 3e8
dz = Lz / Nz
z = range(0, Nz, Nz) .* dz * 1e6

# Calculate number of time steps for one roundtrip
dt = dz / (2 * c0)
#τ = 36     # Choose for testing grid
τ = 512
#τ = 640
# Accounts for reflection from front surface
excess = 0.4e-5 / dz
# Roundtrip time
T = 6 * τ * dt + 2 * ((Nz - excess) * dz / c0)
nsteps = Int(ceil(T / dt))



