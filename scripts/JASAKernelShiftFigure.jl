module JASAKernelShiftFigure

using Plots
using Measures
using Formatting: format

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Bumps.jl")
using .StimGen: hc
using .PeakPick: peakinterp
using .Bumps: gauss

# -------------------- helper functions --------------------
function smooth(xs, ys; σ=8e-5)
	kernel(t, t0, σ) = exp.( - (t.-t0).^2 ./ (2 * σ^2) )
	smoothed = zeros(length(xs))
	for (k, s) in enumerate(xs)
		smoothed[k] = sum(kernel(xs, s, σ) .* ys)
	end
	return smoothed
end


# -------------------- main script --------------------
f0 = 200.0
fs = 100e3
dur = 0.1
num_h = 12 
mth = 3
mt = 1.040
# σ_fac = 0.062
# σ0 = σ_fac / f0
σ0 = 4.5e-4
idx_lo = Int(round(0.004*fs))
idx_hi = Int(round(0.006*fs))

amps = [1.0 for k ∈ 1:num_h]
mistunings = [1.0 for k ∈ 1:num_h]
mistunings[mth] = mt
ϕs = [0.0 for k ∈ 1:num_h]

sac = hc(dur, f0, fs; amps=amps, mistunings=mistunings, phase=ϕs)
t = collect(0:length(sac)-1) ./ fs
smoothed = smooth(t, sac; σ=σ0)
smoothed ./= maximum(smoothed[idx_lo:idx_hi])

# Peak of Smoothed
idx_max = argmax(smoothed[idx_lo:idx_hi]) + idx_lo - 1
peak_smooth = peakinterp(t[idx_max-1:idx_max+1], smoothed[idx_max-1:idx_max+1])
println("Peak of smoothed: $peak_smooth ($(1/peak_smooth) Hz)")

# Peak of sac
idx_max  = argmax(sac[idx_lo:idx_hi]) + idx_lo - 1
peak_sac = peakinterp(t[idx_max-1:idx_max+1], sac[idx_max-1:idx_max+1])
println("Peak of sac: $peak_sac ($(1/peak_sac) Hz)")


## Plot

# mistuned harmonic
p = plot(t, cos.(2π*f0*mth*mt*t) ./ num_h,
		 lw=1.50, lc=:black, la=0.5, label="mistuned component", margin=5mm)

# other harmonics
# harmonics = collect(1:num_h)
# deleteat!(harmonics, mth)
# for k ∈ harmonics 
# 	plot!(p, t, cos.(2π*f0*k*t) ./ num_h, ls=:dash, lc=:black, la=0.5, lw=0.5)
# end

# sac and smoothed
xticks = collect(0.003:0.0005:0.007)
xticks = (xticks, xticks .* 1000)
plot!(p, t, sac,
	  legend=:bottomleft,
	  xlims=(0.00325, 0.00675),
	  ylims=(-0.85, 1.2),
	  lw=2.5,
	  lc=:red,
	  tickfontsize=16,
	  xlabel="Delay (ms)",
	  ylabel="Magnitude (normalized)",
	  guidefontsize=16,
	  size=(700,500),
	  yticks=collect(-0.80:0.2:1.2),
	  xrotation=45,
	  xticks = xticks, 
	  framestyle=:box,
	  label="SACF",
	  legendfontsize=12)
# plot!(p, t, smoothed, lw=2.0, lc=:black, la=0.75, ls=:dash)

# best fitting template
plot!(p, t, gauss(t, peak_smooth, σ0), lw=2.5, lc=1, ls=:dash,
	  label="best template")

# plot peak locations
upoff = 0.05
lwa = 1.75
plot!(p, [peak_smooth, peak_smooth], [0.96, 1.06+upoff], lc=:black, lw=lwa, label="")
plot!(p, [peak_smooth-0.00005, peak_smooth], [1.03+upoff, 1.03+upoff], lc=:black, lw=lwa, label="")
annotate!(p, peak_smooth-0.000090, 1.03+upoff,
		  text( format("{:.3f} (ms)", peak_smooth*1000), 16, :right))
plot!(p, [peak_sac, peak_sac], [0.92, 1.06+upoff], lc=:black, lw=lwa, label="")
plot!(p, [peak_sac+0.00005, peak_sac], [1.03+upoff, 1.03+upoff], lc=:black, lw=lwa, label="")
annotate!(p, peak_sac+0.000090, 1.03+upoff,
		  text( format("{:.3f} (ms)", peak_sac*1000), 16, :left))



end
