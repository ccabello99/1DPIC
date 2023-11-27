using CSV 
using DataFrames
using Plots
using FFTW
include("Grid.jl")
 
# Read the CSVs for incident and reflected fields 
data1 = CSV.read("Ey_inc.csv",DataFrame)
Ey_inc = data1[!, 1]
data2 = CSV.read("Ey_ref.csv",DataFrame)
Ey_ref = data2[!, 1]

# Calculate the number of zeros to pad
zeros_pad1 = nsteps - length(Ey_inc)
zeros_pad2 = nsteps - length(Ey_ref)

# Zero-pad the vector
Ey_inc = vcat(Ey_inc, zeros(zeros_pad1))
Ey_ref = vcat(Ey_ref, zeros(zeros_pad2))

# Time (fs) and frequency (Hz) axis
t = range(1, nsteps, nsteps) .* dt .* 1e15
fs = 1 / dt
f0 = 3e8 / (800e-9)
n = Int(nsteps/2)
f = (fs/nsteps) .* range(0, n, n)

# Spectral Intensity
ft_inc = fft(Ey_inc)
X_inc = abs.(ft_inc)
ft_ref = fft(Ey_ref)
X_ref = abs.(ft_ref)

# Spectral Phase
#ϕ_inc = atan.(imag(ft_inc) / real(ft_inc))
#ϕ_ref = atan.(imag(ft_ref) / real(ft_ref))

# Unwrap phase
#for j in 1:n
#    for k in 1:1000
#       if ϕ_inc[j] < ϕ_inc[j+1]
#            ϕ_inc[j+1] = ϕ_inc[j+1] - π
#        end
#        if ϕ_ref[j] < ϕ_ref[j+1]
#            ϕ_ref[j+1] = ϕ_ref[j+1] - π
#        end
#    end
#end


p = plot(t .- 91.595, 2 .* Ey_inc, xlabel="Time (fs)", ylabel="Electric Field Strength (V/m)", lw = 2, color="blue", label="Incident", title="Incident and Reflected Time-Dependent Fields")
plot!(t .- 64.947, 2 .* Ey_ref, lw = 2, color="red", label="Reflected")
#45.21      # L = 0.1
#55.76      # L = 0.05
#64.947     # L=0.005
#vline!([0], lw=2, color="black")
xlims!(-15, 15)

#p = plot(f ./ f0, [ϕ_inc[1:n] .+ f[1:n] .* .25e-13 .- 0,ϕ_ref[1:n] .+ f[1:n] .* 0.41e-12 .+ 20], xlabel="Frequency (fund. freq.)", ylabel="Spectral Phase (rad)")

p = plot(f ./ f0, X_inc[1:n] ./ maximum(abs.(X_inc)), yaxis=:log, label="Incident", title="Incident and Reflected Field Spectra", xlabel="Frequency (fund. freq.)", ylabel="Spectral Intensity (arb. units)", lw = 2, color="blue")
plot!(f ./ f0, X_ref[1:n] ./ maximum(abs.(X_ref)), yaxis=:log, label="Reflected", lw = 2, color="red")
ylims!(1e-20, 1.5e1)
#xlims!(0, 60)
xlims!(0, 10)
#ylims!(1e-5, 1.5e1)
#xticks!(range(1,20,20))

#p = plot(1e9*c0 ./ f, X_inc[1:n] ./ maximum(abs.(X_inc)), label="Incident Field", title="Incident and Reflected Field Spectra", xlabel="Wavelength (nm)", ylabel="Spectral Intensity (arb. units)")
#plot!(1e9*c0 ./ f, X_ref[1:n] ./ maximum(abs.(X_ref)), label="Reflected Field", title="Incident and Reflected Field Spectra", xlabel="Wavelength (nm)", ylabel="Spectral Intensity (arb. units)")
#ylims=[0, 1]
#xlims!(200, 2000)
#xticks!(range(200,2000,10))

display(p)

# Save vector graphics (pdf)

# For Time Dependent Fields
#savefig(p, "Fields_L0005_a100.pdf")
# For Spectra
savefig(p, "Spect_L0005_a100_first10.pdf")
