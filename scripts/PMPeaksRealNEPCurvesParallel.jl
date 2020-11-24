
## Set up workers
using Distributed

availprocs = length(Sys.cpu_info())
if nprocs() < availprocs
	addprocs(availprocs - nprocs())
end

using Plots
using Suppressor
using AuditoryFilters
using SharedArrays
using Formatting: format
using Statistics: mean, std
using JLD
using Dates

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
@everywhere include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")
include("../util/EdgePitchUtils.jl")
using .EdgePitchUtils: edgepitch

## Load HCC data
hcc_lp = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/hcc_lp_data_clean.jld")
hcc_hp = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/hcc_hp_data_clean.jld")

## Experiment Parameters
fs = 50e3
stim_dur = 0.5  # 0.500
rand_seed = :none
f0s = [ 200.0 * 2 ^ (k/6) for k ∈ -4:22 ] 	# main for figures
# f0s = [ 200.0 * 2 ^ (k/6) for k ∈ -4:4:22 ] 	# main for figures
# f0s = [ 200.0 * 2 ^ (k/3) for k ∈ -1:10 ]
# f0s = [ 200.0 * 2 ^ (k/6) for k ∈ 0:22 ]
# f0s = [ 200.0 * 2 ^ (k/3) for k ∈ -1:5 ]
# f0s = hcc_hp["freqs"]
println(f0s)
num_repeats = 20

num_channels = 50
cf_lo = 100.0
ac_norm = true
ac_dur = stim_dur 
sacparams = PitchModelPeaks.SACParams(num_channels, cf_lo, ac_norm, ac_dur) 
fb = make_erb_filterbank(fs, num_channels, cf_lo)

# σ0 = 4e-4
# skew = 0.03
# τ200 = 0.015
# τ2500 = 0.010
# modelparams = PitchModelPeaks.ModelParams(σ0, τ200, τ2500, skew)
# σ0 = 6.25e-4		#smoothing kernel width at 155.0 Hz
# skew = 0.03 	# Corrective skewing of kernel shape
# τ200 = 0.020	# Integration window lenth at 200 Hz
# τ2500 = 0.010 	# Integration window legnth at 2500 Hz
# bump = :hat
σ0 = 4.50e-4		#smoothing kernel width at 155.0 Hz
skew = 0.030	# Corrective skewing of kernel shape
τ200 = 0.020 	# Integration window lenth at 200 Hz
τ2500 = 0.010 	# Integration window legnth at 2500 Hz
bump = :gauss
modelparams = PitchModelPeaks.ModelParams(σ0, bump, τ200, τ2500, skew)
shthresh = 0.0450	# threshold for declaring subharmonic

## For parallel processing, make closure available
# @everywhere begin 
# 	function find_pitch(stim)
# 		PitchModelPeaks.estimate_τ0(stim, fs; params=modelparams,
# 											  sacparams=sacparams,
# 											  fb=fb)
# 	end
# end

## Plot HCC data
xticks=[collect(100.0:100.0:1000.0); 2000.0; 2500.0]
xticklabels = [format("{:.1f}", tick) for tick ∈ xticks]

p = plot(hcc_lp["freqs"], hcc_lp["mt"], yerr=(hcc_lp["yerrlo"], hcc_lp["yerrhi"]),
		 xscale=:log10, ylims=(0.8, 1.4), lc=:black, la=0.75, markerstrokecolor=:black,
		 xticks=(xticks, xticklabels), label="Humans")
plot!(p, hcc_hp["freqs"], hcc_hp["mt"], yerr=(hcc_hp["yerrlo"], hcc_hp["yerrhi"]),
	  lc=:black, la=0.75, markerstrokecolor=:black, label="Humans")

# Plot HCC Ideal Curves (15, 30, 60 ms)
plot!(p, f0s, (edgepitch.(f0s, 0.015; kind=:lp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5, label=:none)
plot!(p, f0s, (edgepitch.(f0s, 0.015; kind=:hp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5, label=:none)
plot!(p, f0s, (edgepitch.(f0s, 0.030; kind=:lp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5, label=:none)
plot!(p, f0s, (edgepitch.(f0s, 0.030; kind=:hp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5, label=:none)
plot!(p, f0s, (edgepitch.(f0s, 0.060; kind=:lp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5, label=:none)
plot!(p, f0s, (edgepitch.(f0s, 0.060; kind=:hp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5, label=:none)

#kinds = [:lp, :hp]
kinds = [:hp, :lp]
stims_saved = Dict()

for kind ∈ kinds
	## Run trials
	println("-------------------- $kind --------------------")
	ests = []
	mt_μs = []
	mt_stds = []
	mt_μs_hcc = []
	mt_stds_hcc = []
	mt_μs_peak = []
	mt_stds_peak = []
	stds = []
	stims_all = []
	@time for (k, f0) ∈ enumerate(f0s)
		print(format("{:.1f} ($k of $(length(f0s))), of $num_repeats: ", f0))
		stim_len = Int(fs*stim_dur)
		stims = zeros(stim_len, num_repeats) 
		for n ∈ 1:num_repeats
			stims[:,n] = StimGen.nepstim(f0, fs, stim_dur, 0.0; kind=kind, seed=rand_seed)
		end
		# ests_local = pmap(find_pitch, eachcol(stims))
		ests_local = pmap(eachcol(stims)) do stim
			PitchModelPeaks.estimate_τ0(stim, fs; params=modelparams,
												  sacparams=sacparams,
												  fb=fb,
												  returnpeak=true,
												  returnhccest=true,
												  shthresh=shthresh)

		end
		mts = (1 ./ [val[1] for val in ests_local]) ./ f0
		mt_μ = mean(mts)
		mt_std = std(mts)
		push!(mt_μs, mt_μ)
		push!(mt_stds, mt_std)
		
		mts_peak = (1 ./ [val[2] for val in ests_local]) ./ f0
		mt_μ_peak = mean(mts_peak)
		mt_std_peak = std(mts_peak)
		push!(mt_μs_peak, mt_μ_peak)
		push!(mt_stds_peak, mt_std_peak)

		mts_hcc = (1 ./ [val[3] for val in ests_local]) ./ f0
		mt_μ_hcc = mean(mts_hcc)
		mt_std_hcc = std(mts_hcc)
		push!(mt_μs_hcc, mt_μ_hcc)
		push!(mt_stds_hcc, mt_std_hcc)

		push!(ests, mean(1 ./ ests_local[1]))
		push!(stds, std(1 ./ ests_local[1]))

		push!(stims_all, stims)

		println()
	end


	# Plot results over HCC data 
	plot!(p, f0s, mt_μs, yerr=mt_stds, lc=:red, label="model")
	plot!(p, f0s, mt_μs_peak, yerr=mt_stds_peak, lc=:orange, ls=:dash, label="first peak")
	plot!(p, f0s, mt_μs_hcc, yerr=mt_stds_hcc, lc=:blue, ls=:dash, label="hcc (non-smoothed)")

	# Save the stims in environment/module for examination of errors
	stims_saved[kind] = stims_all

	## Save results
	timestamp = Dates.format(Dates.now(), "yy-mm-dd-HH-M")
	JLD.save(
			 "/home/dahlbom/research/AdaptiveSieves/data/latest_$(kind).jld",
			 "f0s", f0s,
			 "mt_μs", mt_μs,
			 "mt_stds", mt_stds,
			 "mt_μs_peak", mt_μs_peak,
			 "mt_stds_peak", mt_stds_peak,
			 "mt_μs_hcc", mt_μs_hcc,
			 "mt_stds_hcc", mt_stds_hcc,
			 "ests", ests,
			 "stds", stds,
			 "τ200", τ200,
			 "τ2500", τ2500,
			 "ac_dur", ac_dur)
	JLD.save(
		 "/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_$(kind)_$timestamp.jld",
		 "f0s", f0s,
		 "mt_μs", mt_μs,
		 "mt_stds", mt_stds,
		 "mt_μs_hcc", mt_μs_hcc,
		 "mt_stds_hcc", mt_stds_hcc,
		 "mt_μs_peak", mt_μs_peak,
		 "mt_stds_peak", mt_stds_peak,
		 "ests", ests,
		 "stds", stds,
		 "τ200", τ200,
		 "τ2500", τ2500,
		 "ac_dur", ac_dur)
end

