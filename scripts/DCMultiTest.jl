module DCMultiTest

using Plots
using AuditoryFilters
using Formatting: format
using Statistics: mean, std
using JLD
using Dates

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")

dc_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/dc_data_clean.jld")

## Experiment Parameters
fs = 50e3
f0 = 155.0
stim_dur = 0.100
mistunings = dc_data["mistunings"]
println(mistunings)
# mistunings = [0.88:0.025:1.12...]
step_size = 0.0025# 0.0025
mt_vals = collect(mistunings[1]:step_size:mistunings[end])
# mt_vals = collect(1.0:step_size:mistunings[end])
num_h = 12
mt_nums = collect(1:5) 
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
# σ0 = 4.0e-4		# Smoothing kernel width at 155.0 Hz
# σ0 = 6.00e-3 	# 4.0e-4 for gauss; 2.00e-3 for cosbump; 3.00e-3 for cosbump2; 1.75e-3 cos0, 0.025
# skew = 1.00 	# Corrective skewing of kernel shape
# τ200 = 0.025 	# Integration window lenth at 200 Hz; 0.02
# τ2500 = 0.014 	# Integration window legnth at 2500 Hz; 0.014
# bump = :cosbump4 # :gauss
# σ0 = 03.00e-4 	# 4.0e-4 for gauss; 2.00e-3 for cosbump; 3.00e-3 for cosbump2; 1.75e-3 cos0, 0.025
# skew = 0.03 	# Corrective skewing of kernel shape
# τ200 = 0.020 	# Integration window lenth at 200 Hz; 0.02
# τ2500 = 0.014 	# Integration window legnth at 2500 Hz; 0.014
# bump = :gauss # :gauss
σ0 = 6.25e-4		#smoothing kernel width at 155.0 Hz
skew = 0.03 	# Corrective skewing of kernel shape
τ200 = 0.020	# Integration window lenth at 200 Hz
τ2500 = 0.010 	# Integration window legnth at 2500 Hz
bump = :hat
modelparams = PitchModelPeaks.ModelParams(σ0, bump, τ200, τ2500, skew)


f0_ests_all = Dict{Int, Array{Float64,1}}()
for mt_num ∈ mt_nums 
	println("-------------------- Partial $mt_num (of $(length(mt_nums))) --------------------")
	f0_ests = Array{Float64,1}() 
	for (k, mt) ∈ enumerate(mt_vals)
		println("$k of $(length(mt_vals))")
		mts = [1.0 for k ∈ 1:num_h]
		mts[mt_num] = mt
		stim = StimGen.hc(stim_dur, f0, fs; amps=amps,
											mistunings=mts,
											phase=ϕs)
		est = PitchModelPeaks.estimate_τ0(stim, fs;
										  params=modelparams,
										  sacparams=sacparams,
										  fb=fb,
										  env=false)
		push!(f0_ests, 1/est)
	end
	f0_ests_all[mt_num] = f0_ests
end


timestamp = Dates.format(Dates.now(), "yy-mm-dd-HH-M")
JLD.save(
		 "/home/dahlbom/research/AdaptiveSieves/data/DCMulti_latest.jld",
		 "f0_ests", f0_ests_all,
		 "mt_vals", mt_vals,
		 "τ200", τ200,
		 "τ2500", τ2500,
		 "ac_dur", ac_dur)
JLD.save(
	 "/home/dahlbom/research/AdaptiveSieves/data/DCMulti_$timestamp.jld",
	 "f0_ests", f0_ests_all,
	 "mt_vals", mt_vals,
	 "τ200", τ200,
	 "τ2500", τ2500,
	 "ac_dur", ac_dur)
end

