module PitchModelPeaks

using AuditoryFilters
using Plots 	# remove after debug

include("../util/StimGen.jl")
using .StimGen: sac
include("../core/Sigmas.jl")
include("../core/Sieves.jl")


function gensac(signal, fs)
	# Params
	num_channels = 50
	lo_freq = 100.0
	norm = true
	dur = 0.100	

	# Action
	fb = make_erb_filterbank(fs, num_channels, lo_freq)	
	sac(signal, fb, fs; dur=dur, norm=norm)			
end

function runmodel(signal, fs)
	# Params
	σ0 = 4e-4
	# sig(t0) = Sigmas.uniform(σ0=σ0, num_bumps=1)
	# sieve(t, t0) = Sieves.sievegauss(t, t0, sig(t0); skew=skew, env=:none)

	# Generate the Summary Autocorrelogram
	sac = gensac(signal, fs)
	τs = collect(0:length(sac)-1) ./ fs

	# Make initial estimate (simple here
	τ_min = (1/3000)
	idcs = findall(x -> x > τ_min, τs)
	idx_max = argmax(sac[idcs])
	τ_est = τs[idcs][idx_max]
	println("initial estimate: ", 1/τ_est)

	# 
end

end
