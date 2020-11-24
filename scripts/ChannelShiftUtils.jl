module ChannelShiftUtils

using FFTW

include("/home/dahlbom/research/AdaptiveSieves/util/ACUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Salience.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Sieves.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Sigmas.jl")

export sum_channels_locally, weight, w_gauss_full, smooth_channels, ac_channels

# -------------------- Functions --------------------
function sum_channels_locally(chans, cfs_in, cfs_out, weightf::Function)
	sig_len, num_chans_in = size(chans)
	chans_summed = zeros(sig_len, length(cfs_out))
	for (n, cf_out) ∈ enumerate(cfs_out)
		for k ∈ 1:num_chans_in
			w = weightf(cf_out, cfs_in[k])
			chans_summed[:,n] .+= w .* chans[:,k]
		end
	end
	chans_summed
end

function weight(cf, f, w)
	f_hi = cf * (2^(w/2))
	f_lo = cf / (2^(w/2))
	# f_hi = cf + 150.0
	# f_lo = cf - 150.0
	if f_lo < f <= f_hi	
		return 1.0
	end
	return 0.0
end

w_gauss_full(cf_sd, cf_an; σ=2) = exp( - abs( log(cf_an) - log(cf_sd) )^2 / σ^2 )

function ac_channels(chans)
	chans_ac = zeros(size(chans))
	chans
	for n ∈ 1:size(chans)[2]
		chans_ac[:,n] = ACUtils.ac_fft(chans[:,n])
	end
	chans_ac
end

function smooth_channels(chans, cfs, sigmaf, fs; bump=:gauss, skew=0.0, env=:none, usefft=false)
	sig_len, num_chans = size(chans)
	ts = collect(0:sig_len-1) ./ fs
	@assert num_chans == length(cfs)
	chans_smoothed = zeros(sig_len, num_chans)
	for n ∈ 1:size(chans)[2]
		σ = sigmaf(cfs[n])
		sig(t0) = Sigmas.uniform(σ0=σ, num_bumps=1)
		sieve(t, t0) = Sieves.sievefunc(t, t0, sig(t0); bump=bump, skew=skew, env=env)
		if usefft
			sieve_vec = sieve.(ts, ts[end]/2)
			norm = sum(abs.(sieve_vec))
			if norm != 0.0
				sieve_vec ./= norm
			end
			chans_smoothed[:,n] = real.(ifft( fft(sieve_vec) .* fft(chans[:,n])))
		else	
			chans_smoothed[:,n] = Salience.salience(ts, ts, chans[:,n], sieve)
		end
	end
	chans_smoothed
end

end
