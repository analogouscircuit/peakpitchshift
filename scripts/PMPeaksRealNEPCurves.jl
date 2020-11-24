module PMPeaksRealNEPCurves

using Plots
using AuditoryFilters
using Formatting: format
using Statistics: mean, std
using JLD
using Dates

include("../core/PitchModelPeaks.jl")
using .PitchModelPeaks: estimate_τ0, SACParams, ModelParams
include("../util/StimGen.jl")


# Experiment Parameters
fs = 50e3
stim_dur = 0.25
f0s = [ 200.0 * 2 ^ (k/6) for k ∈ -4:22 ]
# f0s = [ 200.0 * 2 ^ (k/3) for k ∈ -1:10 ]
println(f0s)
num_repeats = 15

num_channels = 50
cf_lo = 100.0
ac_norm = true
ac_dur = stim_dur 
sacparams = SACParams(num_channels, cf_lo, ac_norm, ac_dur) 
fb = make_erb_filterbank(fs, num_channels, cf_lo)

# σ0 = 4e-4
# τ200 = 0.02
# τ2500 = 0.014
# skew = 0.0
# bump = :gauss # :gauss
# modelparams = ModelParams(σ0, bump, τ200, τ2500, skew)
# σ0 = 6.25e-4		#smoothing kernel width at 155.0 Hz
# skew = 0.03 	# Corrective skewing of kernel shape
# τ200 = 0.020	# Integration window lenth at 200 Hz
# τ2500 = 0.010 	# Integration window legnth at 2500 Hz
# bump = :hat
σ0 = 4.40e-4		#smoothing kernel width at 155.0 Hz
skew = 0.03 	# Corrective skewing of kernel shape
τ200 = 0.020	# Integration window lenth at 200 Hz
τ2500 = 0.010 	# Integration window legnth at 2500 Hz
bump = :gauss
modelparams = PitchModelPeaks.ModelParams(σ0, bump, τ200, τ2500, skew)

kinds = [:hp, :lp]

for kind ∈ kinds
	println("-------------------- $kind --------------------")
	ests = []
	mt_μs = []
	mt_stds = []
	stds = []
	stims_all = []
	for (k, f0) ∈ enumerate(f0s)
		print(format("{:.1f} ($k of $(length(f0s))), of $num_repeats: ", f0))
		ests_local = []
		stims_local = []
		for n ∈ 1:num_repeats
			print("$n ")
			stim = StimGen.nepstim(f0, fs, stim_dur, 0.0; kind=kind)
			est = estimate_τ0(stim, fs;
							  params=modelparams,
							  sacparams=sacparams,
							  fb=fb, 
							  verbose=false,
							  plt=false)
			push!(ests_local, est)
			push!(stims_local, stim)
		end
		mts = (1 ./ ests_local) ./ f0
		mt_μ = mean(mts)
		mt_std = std(mts)
		push!(ests, mean(1 ./ ests_local))
		push!(stds, std(1 ./ ests_local))
		push!(mt_μs, mt_μ)
		push!(mt_stds, mt_std)
		push!(stims_all, stims_local)
		println()
	end

	# p = plot(f0s, ests, yerr=stds, xscale=:log10)
	p = plot(f0s, mt_μs, yerr=mt_stds, xscale=:log10)

	timestamp = Dates.format(Dates.now(), "yy-mm-dd-HH-M")
	save("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_$(kind)_$timestamp.jld",
		 "f0s", f0s,
		 "mt_μs", mt_μs,
		 "mt_stds", mt_stds,
		 "ests", ests,
		 "stds", stds,
		 "τ200", τ200,
		 "τ2500", τ2500,
		 "ac_dur", ac_dur)
end

end
