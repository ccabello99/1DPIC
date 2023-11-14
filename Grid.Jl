using FFTW
# Size of the grid
Nz = 256      # Final number
Lz = 4e-5

# Number of time steps
nsteps = 10000
c0 = 3e8


# Cell size and time stepping
fmax = 1e9
dz = Lz / Nz
dt = dz / (2 * c0)
z = range(0, Nz, Nz) .* dz * 1e6
ks = 2π / dz
kz = fftfreq(Nz, ks)

fmax = 1e17
#nfreq = 2000
#freq = range(-fmax, fmax, length=nfreq)
