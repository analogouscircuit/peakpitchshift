module PMPeaksRealMTPCurves

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

p = plot(dc_data["mistunings"], dc_data["f0s"], yerr=(dc_data["yerrlo"], dc_data["yerrhi"]),
		lc=:black, la=.5)

## Experiment Parameters
fs = 50e3
stim_dur = 0.125
mts = 

# SAC params
num_channels = 50
cf_lo = 100.0
ac_norm = true
ac_dur = stim_dur 
sacparams = PitchModelPeaks.SACParams(num_channels, cf_lo, ac_norm, ac_dur) 
fb = make_erb_filterbank(fs, num_channels, cf_lo)

# Pitch model params
σ0 = 4e-4
τ200 = 0.014
τ2500 = 0.014
skew = 0.0
modelparams = PitchModelPeaks.ModelParams(σ0, τ200, τ2500, skew)


# for (k, mt) ∈ enumerate(mts)
# 
# end


end

