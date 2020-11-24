module JASAMooreDataFig

using Plots 
using JLD
using AuditoryFilters

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")


# Stimuli params
fs = 50e3
f0s = [100, 200, 400]
partials = 1:6
mts = [1.01, 1.02, 1.03, 1.04, 1.06, 1.08]
stim_dur = 0.1
num_h = 12
amps = [1.0 for k ∈ 1:num_h]
phase = [0.0 for k ∈ 1:num_h]

# SAC params
num_channels = 50
cf_lo = 100.0
ac_norm = true
ac_dur = stim_dur 
sacparams = PitchModelPeaks.SACParams(num_channels, cf_lo, ac_norm, ac_dur) 
fb = make_erb_filterbank(fs, num_channels, cf_lo)

# Pitch model params
# σ0 = 4.0e-4		#smoothing kernel width at 155.0 Hz
# skew = 0.03 	# Corrective skewing of kernel shape
# τ200 = 0.020 	# Integration window lenth at 200 Hz
# τ2500 = 0.010 	# Integration window legnth at 2500 Hz
# bump = :gauss
# σ0 = 6.25e-4		#smoothing kernel width at 155.0 Hz
# skew = 0.03 	# Corrective skewing of kernel shape
# τ200 = 0.020	# Integration window lenth at 200 Hz
# τ2500 = 0.010 	# Integration window legnth at 2500 Hz
# bump = :hat
σ0 = 4.50e-4		#smoothing kernel width at 155.0 Hz
skew = 0.030 	# Corrective skewing of kernel shape
τ200 = 0.020 	# Integration window length at 200 Hz
τ2500 = 0.010 	# Integration window length at 2500 Hz
bump = :gauss
modelparams = PitchModelPeaks.ModelParams(σ0, bump, τ200, τ2500, skew)
shthres=0.04



data = Dict()
data_peak = Dict()
stims_all = Dict()
stims = []
count = 0
for f0 ∈ f0s
	global data
	global data_peak
	global stims_all
	global stims
	global count
	data[f0] = Dict()
	data_peak[f0] = Dict()
	stims_all[f0] = Dict()
	for partial in partials
		data[f0][partial] = zeros(length(mts))
		data_peak[f0][partial] = zeros(length(mts))
		stims_all[f0][partial] = [] 
		for (k, mt) ∈ enumerate(mts)
			mistunings = [1.0 for q ∈ 1:num_h]
			mistunings[partial] = mt
			stim = StimGen.hc( stim_dur, f0, fs;
							   amps=amps,
							   mistunings = mistunings,
							   phase=phase )
			push!(stims_all[f0][partial], stim)
			est, est_peak = PitchModelPeaks.estimate_τ0(stim, fs;
											  params=modelparams,
											  sacparams=sacparams,
											  fb=fb,
											  returnpeak=true,
											  shthresh=0.08,
											  env=false)
			data[f0][partial][k] = 100.0*(1/(est*f0) - 1.0)
			data_peak[f0][partial][k] = 100.0*(1/(est_peak*f0) - 1.0)
			push!(stims, stim)
			count += 1
			println("(f0: $f0, partial: $partial, mt: $mt: stim #$count")
		end
	end
end

save("/home/dahlbom/research/AdaptiveSieves/data/paper_data/model_moore_complete.jld",
	 "data_all", data,
	 "data_all_peak", data_peak,
	 "mistunings", mts,
	 "stims", stims_all)


end
