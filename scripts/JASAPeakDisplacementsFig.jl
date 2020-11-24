module JASAPeakDisplacementsFig

using Plots
using Measures
using LaTeXStrings

include("/home/dahlbom/research/AdaptiveSieves/util/EdgePitchUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")

τ200 = 0.0556
τ2500 = 0.0139
τ_max(f) = τ2500 + (τ200 - τ2500) * log(2500.0/f) / log(2500.0/200)

fe = 1000.0
τ_p = 1 / EdgePitchUtils.edgepitch(fe, τ_max(fe); kind=:hp)
fs = 50e3
τs = collect(0.0:(1/fs):0.100)
num_peaks = 20

# Plot HP SAC (sinc) and peak locations
sac = EdgePitchUtils.a_hp(τs, fe)
peaks_sac = PitchModelPeaks.findpeaks(τs, sac, 1/fe, fs, num_peaks)
peaks_ideal = [k/fe for k ∈ 1:num_peaks]
p1 = plot(τs, sac, legend=:none,
		  lw = 2.0,
		  xlims=(0.0, 0.0051),
		  ylims=(0.51, 1.5),
		  guidefontsize=16,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  xlabel=L"$\tau$ (s)",
		  ylabel="SACF Mag.",
		  title="1000 Hz NEP (sinc approx.)",
		  framestyle=:box)
bar_lo = 0.9
bar_hi = 1.3
for k ∈ 1:num_peaks
	plot!(p1, [peaks_sac[k], peaks_sac[k]], [bar_lo, bar_hi], lc=:black, ls=:dash, lw=2.0)
	plot!(p1, [peaks_ideal[k], peaks_ideal[k]], [bar_lo, bar_hi], lc=:red, ls=:dot, lw=2.0)
end
plot!(p1, [τ_p, τ_p], [bar_lo, bar_hi], lc=:black, lw=2.0)


# Plot HP SAC displacements
displacements = 100.0 .* (peaks_sac .- peaks_ideal) .* fe
displacements_scaled = displacements ./ collect(1:length(displacements))
p2 = plot(legend=:none,
		  ylims=(-0.3*100,0.3*100),
		  guidefontsize=16,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  xlabel="Harmonic Number",
		  ylabel="Displacement (%)",
		  framestyle=:box)
for (n, displacement) ∈ enumerate(displacements)
	plot!(p2, [n, n], [0, displacement], lc=:black, line=:stem, mc=:black, marker=:none, lw=3.0)
	plot!(p2, [n+0.25, n+0.25], [0, displacements_scaled[n]], 
		  lc=:red, line=:stem, mc=:red, marker=:none, lw=3.0)
end

# Plot MTHP stimulus and peak locations
fs = 200e3
τs = collect(0.0:1/fs:0.2)
f0 = 200.0
mth = 4
mtp = 1.04
num_h = 12
amps = [1.0 for _ ∈ 1:num_h]
phase = [0.0 for _ ∈ 1:num_h]
mistunings = [1.0 for _ ∈ 1:num_h]
mistunings[mth] = mtp
num_peaks = 20 

sac = StimGen.hc(τs[end], f0, fs; amps=amps, mistunings=mistunings, phase=phase)
peaks_sac = PitchModelPeaks.findpeaks(τs, sac, 1/f0, fs, num_peaks)
peaks_ideal = [k/f0 for k ∈ 1:num_peaks]
p3 = plot(τs, sac,
		  xlims=(0.00490, 0.00510),
		  ylims=(0.81, 1.0),
		  lw=2.0,
		  guidefontsize=16,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  xlabel=L"$\tau$ (s)",
		  ylabel="SACF Mag.",
		  title="200 Hz MTHP",
		  legend=:none,
		  framestyle=:box)

p35 = plot(τs, sac,
		  xlims=(0.00490, 0.00510),
		  ylims=(5.0, 6.0),
		  lw=2.0,
		  guidefontsize=16,
		  legendfontsize=12,
		  label="SACF",
		  grid=false,
		  framestyle=:none)


y_lo = 0.85
y_hi = 0.980
pitch = 1.01*f0
plot!(p3, [peaks_sac[1], peaks_sac[1]], [y_lo, y_hi], lc=:black, lw=2.0, ls=:dash)
plot!(p3, [1/f0, 1/f0], [y_lo, y_hi], lc=:red, ls=:dash, lw=2.0)
plot!(p3, [1/pitch, 1/pitch], [y_lo, y_hi], lc=:black, lw=2.0)
plot!(p35, [peaks_sac[1], peaks_sac[1]], [y_lo, y_hi], lc=:black, lw=2.0, ls=:dash, label="peak")
plot!(p35, [1/f0, 1/f0], [y_lo, y_hi], lc=:red, ls=:dash, lw=2.0, label="reference")
plot!(p35, [1/pitch, 1/pitch], [y_lo, y_hi], lc=:black, lw=2.0, label="pitch")

# Plot MTHP displacements
displacements = 100.0 .* (peaks_sac .- peaks_ideal) .* f0
displacements_scaled = displacements ./ collect(1:length(displacements))
p4 = plot(legend=:none,
		  ylims=(-0.15,0.15),
		  guidefontsize=16,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  xlabel="Harmonic Number",
		  ylabel="Displacement (%)",
		  framestyle=:box)
for (n, displacement) ∈ enumerate(displacements)
	plot!(p4, [n, n], [0, displacement], lw=3.0, lc=:black, line=:stem, mc=:black, marker=:none)
	plot!(p4, [n+0.25, n+0.25], [0, displacements_scaled[n]], 
		  lc=:red, line=:stem, mc=:red, marker=:none, lw=3.0)
end

p45 = plot(legend=true,
		  ylims=(-100,-101),
		  guidefontsize=16,
		  legendfontsize=12,
		  grid=false,
		  framestyle=:none)
for (n, displacement) ∈ enumerate(displacements[1:1])
	plot!(p45, [n, n], [0, displacement],
		  lw=3.0, lc=:black, line=:stem, mc=:black, marker=:none, label="unscaled")
	plot!(p45, [n+0.25, n+0.25], [0, displacements_scaled[n]], 
		  lc=:red, line=:stem, mc=:red, marker=:none, lw=3.0, label="scaled")
end

p = plot(p1, p2, p3, p4, layout=(4,1), size=(700, 900), margin=4mm)

p_nep = plot(p1, p2, layout=(2,1), size=(700, 450), margin=4mm)
p_mth = plot(p3, p4, layout=(2,1), size=(700, 450), margin=4mm)

l = @layout [a{0.425w} b{0.425w} c; d{0.425w} d{0.425w} f]
p_horiz = plot(p1, p3, p35, p2, p4, p45,
			   layout=l,
			   size=(1400, 450),
			   margin=5mm)

end
