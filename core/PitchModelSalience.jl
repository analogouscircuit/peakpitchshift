module PitchModelSalience


# -------------------- Packages and Modules --------------------

using AuditoryFilters

include("../util/StimGen.jl")
include("../core/Sigmas.jl")
include("../core/Sieves.jl")
include("../core/Envelopes.jl")
include("../core/Salience.jl")
include("../util/EdgePitchUtils.jl")
include("../util/PeakPick.jl")



# -------------------- Types --------------------

mutable struct SACParams
	num_channels::Int16 	# number of auditory channels
	lo_freq::Float64  		# low cut-off frequency of GT filters
	norm::Bool 				# normalize the SAC
	dur::Float64 			# length of SAC window
end

SACParams() = SACParams(50, 100.0, true, 0.100) 


mutable struct ModelParams
	σ0::Float64 			# width of smoothing kernel at 155 Hz
	τ200::Float64 			# max integration window at 200 Hz
	τ2500::Float64 			# max integration window at 2.5 kHz
	skew::Float64 			# skew factor (distortion of smoothing kernel)
	τ_min::Float64 			# don't consider pitches above 3 kHz
end

ModelParams() = ModelParams(4e-4, 0.05, 0.05, 0.0, 1/2500)



# -------------------- Functions --------------------

"""
	gensac(signal, fs; params=SACParams())

Generate a summary autocorrelagram from an input signal.
"""
function gensac(signal, fs; params=SACParams())
	fb = make_erb_filterbank(fs, params.num_channels, params.lo_freq)	
	StimGen.sac(signal, fb, fs; dur=params.dur, norm=params.norm)			
end


"""
	gensac(signal, fs; params=SACParams())

Generate a summary autocorrelagram from an input signal. Takes a filterbank (do avoid repeated calculation.
"""
function gensac(signal, fs, fb; params=SACParams())
	StimGen.sac(signal, fb, fs; dur=params.dur, norm=params.norm)			
end


"""
	 model_peak(signal, fs)

Run the pitch model (peak picking variant)
"""
function estimate_τ0(signal, fs; params=ModelParams(), verbose=false)
	# Generate the Summary Autocorrelogram
	sac = gensac(signal, fs)
	τs = collect(0:length(sac)-1) ./ fs
	
	# Make initial estimate (simple peak picking here)
	idx_first = searchsortedfirst(τs, params.τ_min)
	idx_max = argmax(sac[idx_first:end])
	println(1/τs[idx_first])
	τ_est = τs[idx_first:end][idx_max]
	
	# Set up helper functions
	maxpeaks = Int(floor(Envelopes.hcc_maxdelay(τ_est; 
												τ200=params.τ200,
												τ2500=params.τ2500) / τ_est))
	sig(t0) = Sigmas.uniform(σ0=params.σ0, num_bumps=maxpeaks)
	env(t, t0) = Envelopes.gamma(t, t0; peakstart = false)
	sieve(t, t0) = Sieves.sievegauss(t, t0, sig(t0); skew=params.skew, env=env)
	peakf(t, t0, x) = Salience.peak_sal(t, t0, x; sievef=sieve)

	# Calculate pitch estimate
	τ0 = Salience.peak_sal(τs, τ_est, sac; sievef=sieve)

	if verbose
		println("Max peaks: $maxpeaks")
		println("Initial guess: $τ_est s ($(1/τ_est) Hz)")
		println("Final estimate: $τ0 s ($(1/τ0) Hz)")
	end

	return τ0
end


end
