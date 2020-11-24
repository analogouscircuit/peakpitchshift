module PACTuningCurve

# -------------------- Packages and Modules --------------------
using AuditoryFilters
using DSP
using Plots
using Plots.Measures
using LinearAlgebra: eachcol

include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")
using .PeakPick: findpeaknear, peakinterp
include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
using .StimGen: hc
include("/home/dahlbom/research/AdaptiveSieves/util/ACUtils.jl")
using .ACUtils: ac_true, ac_fft


# -------------------- Functions --------------------
function peakest(x, y, guess)
	i = findpeaknear(x, y, guess)
	(i == 1) && (i = 2)
	(i == length(y)) && (i = length(y) - 1)
	peakinterp(x[i-1:i+1], y[i-1:i+1])
end

function pacgen(channels; thresh = 0.1)
	e_vals = [mapreduce(x -> x^2, (+), stim) for stim ∈ channels]
	e_max = maximum(e_vals)
	idcs = findall( x -> x > thresh*e_max, e_vals)
	pac_channels = channels[idcs]
	pac = pac_channels[1]
	for k ∈ 2:length(pac_channels)
		pac .*= pac_channels[k]
	end
	return pac
end

function hwr(sig)
	map( x -> x < 0.0 ? 0.0 : x, sig )
end


# -------------------- Trial Parameters --------------------
fs = 100e3
f0 = 155.0 
dur = 0.1
num_chan = 50
num_h = 12
mistuned_h = 4
amps = [1.0 for k ∈ 1:num_h]
ϕs = [0.0 for k ∈ 1:num_h]
mt_vals = collect(0.9:0.01:1.1)
fb = make_erb_filterbank(fs, num_chan, 50.0)


# -------------------- Run Simulation --------------------
ests = []
pacs = []
for mt ∈ mt_vals
	println("$mt")
	mistunings = [1.0 for k ∈ 1:num_h]
	mistunings[mistuned_h] = mt
	stim = hc( dur, f0, fs; amps=amps, mistunings=mistunings, phase=ϕs)
	t = collect(0:length(stim)-1) ./ fs
	stim_filt = filt(fb, stim)
	stim_chans = [stim for stim ∈ eachcol(stim_filt)]
	stim_hwr = [hwr(stim) for stim ∈ stim_chans]
	stim_ac = [ac_fft(stim; norm=false) for stim ∈ stim_hwr]
	pac = pacgen(stim_ac)
	push!(pacs, pac)
	peak = peakest(t, pac, 1.0/f0)
	push!(ests, 1/peak)
end


# -------------------- Plot Results --------------------
p = plot(mt_vals, ests)

end
