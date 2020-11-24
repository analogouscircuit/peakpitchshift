module JASAPitchDemoFigure

using Plots
using AuditoryFilters

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")
include("../util/EdgePitchUtils.jl")
using .EdgePitchUtils: edgepitch


## Stimulus parameters
fs = 50e3
stim_dur = 0.1
f0 = 200.0
kind=:hp

## Model Parameters
num_channels = 50
cf_lo = 100.0
ac_norm = true
ac_dur = stim_dur 
sacparams = PitchModelPeaks.SACParams(num_channels, cf_lo, ac_norm, ac_dur) 
fb = make_erb_filterbank(fs, num_channels, cf_lo)

σ0 = 4.5e-4
skew = 0.03
τ200 = 0.020
τ2500 = 0.010
bump = :gauss
modelparams = PitchModelPeaks.ModelParams(σ0, bump, τ200, τ2500, skew)


## Generate Stimulus
# stim = StimGen.nepstim(f0, fs, stim_dur, 0.0; kind=kind)
mistunings = [1.0 for k ∈ 1:12]
mistunings[3] = 1.04
stim = StimGen.hc(stim_dur, f0, fs; amps=[1.0 for k ∈ 1:12],
									 mistunings=mistunings,
									 phase=[0.0 for k ∈ 1:12])


## Run the pitch model
τ, p = PitchModelPeaks.estimate_τ0(stim, fs;
								   params=modelparams,
								   sacparams=sacparams,
								   fb=fb,
								   plt=true)


end
