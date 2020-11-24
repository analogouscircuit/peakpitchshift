module PMPeaksRealDCCurves

using Plots
using AuditoryFilters
using Formatting: format
using Statistics: mean, std
using JLD
using Dates

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")

## Load Darwin-Ciocca Data
dc_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/dc_data_clean.jld")

xaxis = ("Mistuning Factor of 4h Harmonic", font(16))
yaxis = ("Perceived F0 (Hz)", font(16))

p = scatter(dc_data["mistunings"], dc_data["f0s"],
            color=:gray,
            yerror = (dc_data["yerrlo"], dc_data["yerrhi"]),
            markershape=:square,
            markersize=8,
            markerstrokewidth=0.5,
            markerstrokecolor=:black,
            markeralpha=0.3,
            markercolor=:black,
			legend=:none,
			size=(800, 600),
			xaxis=xaxis,
			yaxis=yaxis,
			xtickfontsize=12,
			ytickfontsize=12)

plot!(p, dc_data["mistunings"], dc_data["f0s"], yerr=(dc_data["yerrlo"], dc_data["yerrhi"]),
	 lc=:black, la=.3, ls=:dash, legend=:none)

## Experiment Parameters
fs = 50e3
f0 = 155.0
stim_dur = 0.100
mistunings = dc_data["mistunings"]
step_size = 0.0025
mt_vals = collect(mistunings[1]:step_size:mistunings[end]+0.01)
num_h = 12
mt_num = 4
amps = [1.0 for k ∈ 1:num_h]
ϕs = [-π/2 for k ∈ 1:num_h]

# SAC params
num_channels = 50
cf_lo = 100.0
ac_norm = true
ac_dur = stim_dur 
sacparams = PitchModelPeaks.SACParams(num_channels, cf_lo, ac_norm, ac_dur) 
fb = make_erb_filterbank(fs, num_channels, cf_lo)

# Pitch model params
σ0 = 4.50e-4		#smoothing kernel width at 155.0 Hz
skew = 0.030	# Corrective skewing of kernel shape
τ200 = 0.020 	# Integration window lenth at 200 Hz
τ2500 = 0.010 	# Integration window legnth at 2500 Hz
bump = :gauss
# σ0 = 5.0e-4		#smoothing kernel width at 155.0 Hz
# skew = 0.03 	# Corrective skewing of kernel shape
# τ200 = 0.015 	# Integration window lenth at 200 Hz
# τ2500 = 0.010 	# Integration window legnth at 2500 Hz
# bump = :gauss
# modelparams = PitchModelPeaks.ModelParams(σ0, τ200, τ2500, skew)
# σ0 = 6.25e-4		#smoothing kernel width at 155.0 Hz
# skew = 0.30 	# Corrective skewing of kernel shape
# τ200 = 0.025	# Integration window lenth at 200 Hz
# τ2500 = 0.010 	# Integration window legnth at 2500 Hz
# bump = :hat
modelparams = PitchModelPeaks.ModelParams(σ0, bump, τ200, τ2500, skew)
shthresh = 0.04


f0_ests = []
f0_hcc = []
f0_peaks = []
stims = []
for (k, mt) ∈ enumerate(mt_vals)
	println("$k of $(length(mt_vals))")
	mts = [1.0 for k ∈ 1:num_h]
	mts[mt_num] = mt
	stim = StimGen.hc(stim_dur, f0, fs; amps=amps,
					  					mistunings=mts,
										phase=ϕs)
	push!(stims, stim)
	est, peak, est_hcc = PitchModelPeaks.estimate_τ0(stim, fs;
													  params=modelparams,
													  sacparams=sacparams,
													  fb=fb,
													  returnpeak=true,
													  returnhccest=true,
													  env=false,
													  shthresh=shthresh)
	push!(f0_ests, 1/est)
	push!(f0_hcc, 1/est_hcc)
	push!(f0_peaks, 1/peak)
end

plot!(p, mt_vals, f0_ests, lc=:red, legend=:none, lw=2.0)
plot!(p, mt_vals, f0_peaks, lc=:blue, legend=:none, ls=:dot, lw=1.5)
plot!(p, mt_vals, f0_hcc, lc=:black, legend=:none, ls=:dashdot, lw=1.5)


timestamp = Dates.format(Dates.now(), "yy-mm-dd-HH-M")
JLD.save(
		 "/home/dahlbom/research/AdaptiveSieves/data/DC_latest.jld",
		 "f0s_ests", f0_ests,
		 "mt_vals", mt_vals,
		 "f0_peaks", f0_peaks,
		 "f0_hcc", f0_hcc,
		 "τ200", τ200,
		 "τ2500", τ2500,
		 "ac_dur", ac_dur)
JLD.save(
	 "/home/dahlbom/research/AdaptiveSieves/data/DC_$timestamp.jld",
	 "f0s_ests", f0_ests,
	 "mt_vals", mt_vals,
	 "f0_peaks", f0_peaks,
	 "f0_hcc", f0_hcc,
	 "τ200", τ200,
	 "τ2500", τ2500,
	 "ac_dur", ac_dur)
end

