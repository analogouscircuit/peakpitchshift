module PitchModelPeaks

export estimate_τ0, SACParams, ModelParams

# -------------------- Packages and Modules --------------------

using AuditoryFilters
using Distributed
using Suppressor
using Plots
using Measures
using Formatting: format

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Sigmas.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Sieves.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Envelopes.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Salience.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/EdgePitchUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")


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
	bump::Symbol 			# :gauss, :cosbump, :cosbump2, :rect, :laplace, :hat, :powgauss
	τ200::Float64 			# max integration window at 200 Hz
	τ2500::Float64 			# max integration window at 2.5 kHz
	skew::Float64 			# skew factor (distortion of smoothing kernel)
end

ModelParams() = ModelParams(4e-4, :gauss, 0.025, 0.010, 0.0)


# -------------------- Utility Functions --------------------

"""
	findpeak(τs, sac)

Finds all pos->neg zero crossings of gradient.  Selects the maximum of these as
the peak.  If no such maxima available, returns maximum of range.
"""
function findpeak(x, y; idx=false)
	@assert length(x) == length(y) "x and y not the same length!"
	g = PeakPick.grad1D(y, x[2] - x[1])
	cross_idcs = []
	for k ∈ 2:length(y)-1
		if g[k-1] >= 0 && g[k+1] < 0
			push!(cross_idcs, k)
		end
	end
	if length(cross_idcs) == 0
		idx_max = argmax(y)
		peak = x[ idx_max ]
		truemax = false
	else
		idx_max = cross_idcs[ argmax(y[cross_idcs]) ]
		peak = PeakPick.peakinterp(x[idx_max-1:idx_max+1], y[idx_max-1:idx_max+1])
		truemax = true
	end
	if idx == true
		return peak, truemax, idx_max
	end
	peak, truemax
end

"""
	findpeaks(τs, smoothed, τ_est, fs, maxpeaks)

Find peaks for smoothed SAC. Must be given an initial estimate, then picks 
subsequent peaks sequentially based on τ0 associated with last peak picked.
"""
function findpeaks(τs, smoothed, τ_est, fs, maxpeaks)
	half_width_t = τ_est/2
	half_width_n = Int(floor(half_width_t * fs ))
	peaks = zeros(maxpeaks)
	i_c = Int(floor(τ_est*fs)) 	# center index
	#i_lo = i_c - Int(floor(half_width_n/2))
	i_lo = i_c - half_width_n 
	i_lo = i_lo < 1 ? 1 : i_lo
	i_hi = i_c + half_width_n
	i_hi = i_hi > length(smoothed) ? length(smoothed) : i_hi
	peaks[1], truemax = findpeak(τs[i_lo:i_hi], smoothed[i_lo:i_hi])

	# Find subsequent peaks
	for n ∈ 2:maxpeaks
		τ_jump = sum( peaks[1:n-1] .* [1/q for q ∈ 1:(n-1)] )  ./ (n-1) 
		# println("τ_jump avg method: $τ_jump")
		# τ_jump = peaks[n-1] / (n-1)
		# println("τ_jump: $τ_jump")
		half_width_n = Int(floor(fs * τ_jump/4))
		# println("half width steps: $half_width_n")
		τ_c = τ_jump + peaks[n-1]
		# println("τ_c: $τ_c") 
		i_c = Int(floor(τ_c*fs)) 	# center index
		# println("i_c: $i_c")
		i_lo = i_c - half_width_n
		i_lo = i_lo < 1 ? 1 : i_lo
		i_hi = i_c + half_width_n
		i_hi = i_hi > length(smoothed) ? length(smoothed) : i_hi
		# println("idcs: $i_lo, $i_hi")
		# i_peak = argmax(smoothed[i_lo:i_hi]) + i_lo - 1
		# println("i_peak: $i_peak")
		# peaks[n] = PeakPick.peakinterp(τs[i_peak-1:i_peak+1],
		# 							   smoothed[i_peak-1:i_peak+1], testpeak=true)
		# peaks[n] = τs[i_peak]
		peaks[n], truemax = findpeak(τs[i_lo:i_hi], smoothed[i_lo:i_hi])
		if peaks[n] < peaks[n-1]
			println("you messed up picking your first peak, I think.")
			peaks[n] = 2*peaks[n-1]
		end
		# println("Peak n: $(peaks[n])")
	end
	return peaks
end

"""
	findpeaksgrad(τs, smoothed, τ_est, fs, maxpeaks)

Find peaks for smoothed SAC. Must be given an initial estimate, then picks 
subsequent peaks sequentially based on τ0 associated with last peak picked.
"""
function findpeaksgrad(τs, smoothed, τ_est, fs, maxpeaks)
	cross_idcs = PeakPick.crossings(PeakPick.grad1D(smoothed, τs[2] - τs[1]))
	peaks = τs[cross_idcs][2:maxpeaks+1]
	return peaks
end

"""
    firsttrough(sac, fs)

Finds the index of the first trough in a SAC.  Smooths, takes gradient and
finds first zero crossing.
"""
function firsttrough(τs, sac; σ=8e-5, returnsmoothed=false)
	idx_last = searchsortedfirst(τs, 1/50)
	smoothed = smooth(τs[1:idx_last], sac[1:idx_last]; σ=σ)	
	gradient = PeakPick.grad1D(smoothed, τs[2] - τs[1])
	cross_idx = 0
	for idx ∈ 2:length(gradient)-1 
		if gradient[idx-1] <= 0 && gradient[idx+1] > 0
			cross_idx = idx
			break
		end
	end
	if returnsmoothed
		return cross_idx, smoothed
	end
	return cross_idx
end


"""
    firsttrough(sac, fs)

Finds the index of the first trough in a SAC.  Smooths, takes gradient and
finds first zero crossing.
"""
function firsttroughsm(τs, sac; σ=8e-5)
	idx_last = searchsortedfirst(τs, 1/50)
	smoothed = smooth(τs, sac; σ=σ)	
	gradient = PeakPick.grad1D(smoothed[1:idx_last], τs[2] - τs[1])
	cross_idx = 0
	for idx ∈ 2:length(gradient)-1 
		if gradient[idx-1] <= 0 && gradient[idx+1] > 0
			cross_idx = idx
			break
		end
	end
	return cross_idx, smoothed
end

"""
    firstpeak(sac, fs)

Finds the index of the first peak in a SAC.  Smooths, takes gradient and
finds first zero crossing.
"""
function firstpeak(τs, sac; σ=8e-5)
	smoothed = smooth(τs, sac)	
	gradient = PeakPick.grad1D(smoothed, τs[2] - τs[1])
	cross_idx = 0
	for idx ∈ 2:length(gradient)-1 
		if gradient[idx-1] >= 0 && gradient[idx+1] < 0
			cross_idx = idx
			break
		end
	end
	return cross_idx
end

# """
# 	findfirstpeak(sac, fs)
# """
# function findfirstpeak(τs, sac; σ=8e-5)
# 	smoothed = smooth(τs, sac)	
# 	gradient = PeakPick.grad1D(smoothed, τs[2] - τs[1])
# 	# find first trough
# 	cross_idx1 = 0
# 	for idx ∈ 2:length(gradient)-1 
# 		if gradient[idx-1] <= 0 && gradient[idx+1] > 0
# 			cross_idx1 = idx
# 			break
# 		end
# 	end
# 	# find subsequent peak
# 	cross_idx2 = cross_idx1
# 	for idx ∈ cross_idx1:length(gradient)-1 
# 		if gradient[idx-1] >= 0 && gradient[idx+1] < 0
# 			cross_idx2 = idx
# 			break
# 		end
# 	end
# 	return cross_idx2
# end

"""
	smooth(x, y; σ=8e-5)

Smooth with a gaussian kernel. Assumes uniform sampling in xs
"""
function smooth(xs, ys; σ=5-e5) 
	dx = xs[2] - xs[1]
	kernel(t, t0, σ) = exp.( - (t.-t0).^2 ./ (2 * σ^2) )
	smoothed = zeros(length(xs))
	idx_offset = Int(round(2σ/dx))
	for (k, x) in enumerate(xs)
		i_c = Int(round((x - xs[1])/dx)) + 1
		i_lo = i_c - idx_offset 
		(i_lo < 1) && (i_lo = 1)
		i_hi = i_c + idx_offset 
		(i_hi > length(ys)) && (i_hi = length(ys))
		smoothed[k] = sum(kernel(xs[i_lo:i_hi], x, σ) .* ys[i_lo:i_hi])
	end
	return smoothed
end


"""
	gensac(signal, fs; params=SACParams())

Generate a summary autocorrelagram from an input signal. Generates own
filterbank if not provided. Can provide filter bank if going to run multiple
times and want to avoid overhead of recalculating filters.
"""
function gensac(signal, fs; params=SACParams(), fb=:none, channels=false, acchannels=false)
	if fb == :none
		fb = make_erb_filterbank(fs, params.num_channels, params.lo_freq)	
	end
	StimGen.sac(signal, fb, fs; 
			    dur=params.dur,
			    norm=params.norm,
			    channels=channels,
			    acchannels=acchannels)			
end


# -------------------- Model --------------------

"""
	estimate_τ0(signal, fs; params=ModelParams(),
							sacparams = SACParams(),
							fb=:none,
							verbose=false, 
							plt=false)

Run the pitch model (peak picking variant). Generates a summary autocorrelogram,
(SAC); makes an initial pitch estimate; smooths the SAC based on initial
estimate; picks peaks, the number of which is determined by the HCC formula; 
finally estimates pitch based on minimizing L2 distance from ideal harmonic
sieves. Optional plotting to see output of each stage.
"""
function estimate_τ0(signal, fs; params=ModelParams(),
					 			 sacparams = SACParams(),
								 fb=:none,
								 env=false,
								 verbose=false, 
								 plt=false,
								 returnpeak=false,
								 returnhccest=false,
								 shthresh=0.05)

	## Generate the Summary Autocorrelogram (SAC)
	@assert length(signal)/fs >= sacparams.dur
		"Not enough signal to make an autocorrelation of specified length!"
	sac = gensac(signal, fs;
			     params=sacparams,
			     fb=fb,
			     channels=false,
			     acchannels=false)
	# sac = gensac(signal, fs; params=sacparams, fb=fb)
	# sac_squashed = copy(sac)
	τs = collect(0:length(sac)-1) ./ fs

	## Make initial estimate (maximum amp post-trough down to 50 Hz)
	#idx_first, sm_init = firsttrough(τs, sac; returnsmoothed=true)
	idx_first, sm_init = firsttroughsm(τs, sac)
	verbose && println("Not less than: $(τs[idx_first]) ($(1/τs[idx_first]))")
	offset_factor_1 = 1.00 	# in case of early trough
	idx_start = Int(round(offset_factor_1*idx_first))
	# offset_factor_2 = 7 	# shouldn't catch more than one subharmonic
	# idx_stop = Int(round(offset_factor_2*idx_first))
	τ_max = 1.0/80 # lowest pitch to consider
	i_max_delay = Int(round(τ_max*fs))
	idx_stop = i_max_delay
	# sac[1:idx_start-1] .= 0.0 	# kill zero-delay bump
	# squash = ones(length(sac)) .* 0.5
	# squash[1:idx_stop] .= collect(range(1.0, 0.75, length=idx_stop))
	#sac_squashed = sac .* squash 
	sac_squashed = copy(sac)
	idx_max = argmax(sac_squashed[idx_start:idx_stop]) + idx_start - 1 
	verbose && println("First guess: $(τs[idx_max]) ($(1/τs[idx_max]))")
	τ_init = τs[idx_max]

	# Check if it's a subharmonic
	for sub ∈ [3,2]
		idx_half_max = Int(round(idx_max/sub))
		if idx_half_max > idx_start
			half_width = Int(round(idx_half_max/3))
			idx_lo = idx_half_max-half_width # 30
			(idx_lo < 1) && (idx_lo = idx_half_max)
			idx_hi = idx_half_max+half_width # 30
			(idx_hi > length(sac)) && (idx_hi = length(sac))
			# peak_cand, truemax, idx_cand = findpeak(τs[idx_lo:idx_hi],
			# 										sac_squashed[idx_lo:idx_hi]; idx=true)
			peak_cand, truemax, idx_cand = findpeak(τs[idx_lo:idx_hi],
													sm_init[idx_lo:idx_hi]; idx=true)
			idx_cand += idx_lo - 1
			#if (abs(sac_squashed[idx_cand]/sac_squashed[idx_max] - 1) < shthresh
			if (abs(sm_init[idx_cand]/sm_init[idx_max] - 1) < shthresh
				&& truemax
				&& idx_cand > idx_first)
				println("chose a subharmonic at $(τs[idx_max]), ($(1/τs[idx_max]))...")
				idx_max = idx_cand
				break
			end
		end
	end

	# make final estimate
	sac[1:idx_start-1] .= 0.0 	# kill zero-delay bump
	τ_est = PeakPick.peakinterp(τs[idx_max-1:idx_max+1], sac[idx_max-1:idx_max+1], testpeak=true)

	## Set up helper functions (smoothing kernel apparatus)
	#σ0 = params.σ0*155.0*τ_est 	# determine smoothing kernel width
	#σ0 = params.σ0*(1.0 + log(155.0*τ_est)) 	# determine smoothing kernel width
	σ0 = params.σ0*sqrt(155.0*τ_est) 	# determine smoothing kernel width
	# σ0 = params.σ0
	sig(t0) = Sigmas.uniform(σ0=σ0, num_bumps=1)
	# sieve(t, t0) = Sieves.sieve(t, t0, sig(t0); bump=params.bump, skew=params.skew, env=:none)
	sieve(t, t0) = Sieves.sievefunc(t, t0, sig(t0); bump=params.bump, skew=params.skew, env=:none)

	## Get smoothed SAC
	smoothed = Salience.salience(τs, τs, sac, sieve)

	## Find the peaks (maxpeaks of them)
	maxpeaks = Int(floor(Envelopes.hcc_maxdelay(τ_est;
												τ200=params.τ200,
												τ2500=params.τ2500) / τ_est))
	# println("About to find peaks...")
	peaks = findpeaks(τs, smoothed, τ_est, fs, maxpeaks)
	# τ_est = peaks[1] # first peak only from smoothed -- remove this later
	if env
		#envf(t) = Envelopes.gamma(t, τ_est)
		τ_decay = maxpeaks * τ_est
		envf(t) = Envelopes.expenv(t, τ_decay*3.0)
		# envf(t) = Envelopes.tempift(t, τ_est; num_h=5, α=0.031)
		τ0 = EdgePitchUtils.delayfrompeaks(peaks; maxpeaks=Int(maxpeaks*2), envf=envf)
	else
		τ0 = EdgePitchUtils.delayfrompeaks(peaks; maxpeaks=maxpeaks)
	end

	if returnhccest
		# returns same calculation but without using a smoothed SAC
		peaks_hcc = findpeaks(τs, sac, τ_est, fs, maxpeaks)
		τ0_hcc = EdgePitchUtils.delayfrompeaks(peaks_hcc; maxpeaks=maxpeaks)
	end

	## Report and return
	if verbose
		println("σ0: $σ0")
		println("Max peaks: $maxpeaks")
		println("Initial guess: $τ_est s ($(1/τ_est) Hz)")
		# println("Initial guess index: $(idx_first + idx_peak)")
		println("Final estimate: $τ0 s ($(1/τ0) Hz)")
	end

	if plt
		xticks = 0:0.01:0.03
		xticks = (xticks, 1000.0 .* xticks)
		p1 = plot(collect(0:length(signal)-1) ./ fs, signal, 
				  legend=:none,
				  xlabel="time (s)",
				  ylabel="magnitude",
				  title="Input Signal")
		p2 = plot(τs, sac_squashed, 
				  legend=:none,
				  framestyle=:box,
				  lw=2.0,
				  xlabel="delay (ms)",
				  ylabel="magnitude",
				  tickfontsize=14,
				  guidefontsize=16,
				  title="Summary Autocorrelagram Function (SACF)",
				  xticks=xticks,
				  xlims=(0,0.03))
		plot!(p2, [τ_est, τ_est], [0,1], lc=:black, ls=:dash, lw=1.5, la=0.5)
		annotate!(p2, τ_est+0.00025, 0.90,
				  text(format("{:.3f}\n(ms)", τ_est*1000), 12, :italic, :left))
		annotate!(p2, 0.0210, 0.6, text(format("Pitch: {:.2f} Hz", 1/τ_est), 14, :left, :red))
		p25 = plot(τs[1:length(sm_init)], sm_init ./ maximum(sm_init),
				   legend=:none,
				  framestyle=:box,
				   xlabel="delay (s)",
				   ylabel="magnitude",
				   title="Smoothed for initial peak",
				   xlims=(0,0.02))
		plot!(p25, [τ_est, τ_est], [0,1])
		p3 = plot(τs, smoothed,
				  legend=:none,
				  framestyle=:box,
				  xlabel="delay (ms)",
				  ylabel="magnitude",
				  tickfontsize=14,
				  guidefontsize=16,
				  title="Smoothed SACF",
				  lw=2.0,
				  #xlims=(0, min(peaks[end]*2, sacparams.dur)))
				  # xlims=(0, sacparams.dur))
				  xticks=xticks,
				  xlims=(0,0.03))
		for (k, peak) ∈ enumerate(peaks)
			plot!(p3, [peak, peak], [0, 1], lc=:black, ls=:dash, lw=1.5, la=0.5)
			annotate!(p3, peak+0.00025, 0.90, text(format("{:.3f}\n(ms)", peak*1000), 12, :left))
			# plot!(p3, [k*τ_est, k*τ_est], [0, 1], lc=:green, la=0.5)
		end
		annotate!(p3, 0.0210, 0.6, text(format("Pitch: {:.2f} Hz", 1/τ0), 14, :left, :red))
		# plot!(p3, [τ0, τ0], [0, 1], lc=:red)
		# annotate!(p3, τ0, 1.2,
		# 		  text(format("τ0: {:.4f} ms,  f0: {:.2f} Hz", τ0*1000, 1/τ0),
		# 			   12, :left))
		# p = plot(p1, p2, p25, p3, layout=(4,1), legend=:none)
		p = plot(p2, p3, layout=(2,1),
				 legend=:none,
				 size=(700,650),
				 margin=5mm)

		return τ0, p
	end

	if returnpeak && returnhccest
		return τ0, τ_est, τ0_hcc
	elseif returnhccest	
		return τ0, τ0_hcc
	elseif returnpeak
		return τ0, τ_est
	else
		return τ0
	end
end


end
