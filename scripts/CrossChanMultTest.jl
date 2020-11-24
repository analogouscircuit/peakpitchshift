module CrossChanMultTest

using AuditoryFilters
using Plots

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/ACUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/EdgePitchUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")

fs = 50e3
f0 = 200.0
dur = 0.15
fb = make_erb_filterbank(fs, 50, 50.0)
cfs_an = fb.center_frequencies
weighted = true

## Generate signal and HWR-ed filter channels
mistunings = [1.0 for k ∈ 1:12]
mistunings[3] = 1.00
stim = StimGen.hc(dur, f0, fs; mistunings=mistunings)
stim_filt = filt(fb, stim)
stim_hwr = map(y->y < 0.0 ? 0.0 : y, stim_filt)
# stim_hwr = map(y-> y+0.1, stim_hwr)
# sacf = sum(stim_hwr, dims=2)[:,1]

## Weighting function
σ = 1.0 
w(f_an, f_sd) = exp( - abs(log(f_an) - log(f_sd))^2 / (σ^2))
cfs_sd = collect( exp.( range(log(50.0), log(3000.0), length=50) ) )

if weighted
	sd_channels = ones(size(stim_hwr)[1], length(cfs_sd))
	for (n, cf_sd) ∈ enumerate(cfs_sd)
		for (p, cf_an) ∈ enumerate(cfs_an)
			weight = w(cf_an, cf_sd)
			if weight >= 0.1
				sd_channels[:,n] = (stim_hwr[:,p] .* sd_channels[:,n]) *  weight
			end
		end
	end

	## Calculate autocorrelation of each channel
	sd_channels_ac = ones(size(stim_hwr)[1], length(cfs_sd))
	for (n, chan) ∈ enumerate(eachcol(sd_channels))
		sd_channels_ac[:,n] = ACUtils.ac_fft(chan)
	end
	## Sum results
	spacf = sum(sd_channels_ac, dims=2)[:,1]
else
	## All channels multiplied
	spacf = reduce( (x, y) -> x .* y, eachcol(stim_hwr))
	spacf = ACUtils.ac_fft(spacf)
end

# Analysis and plotting
spacf ./= maximum(spacf)
ts = collect(0:length(spacf)-1) ./ fs
peaks = PitchModelPeaks.findpeaks(ts, spacf, 1/f0, fs, 5)
println(1 ./ peaks)
τ0 = EdgePitchUtils.delayfrompeaks(peaks; maxpeaks=3)
println(1/τ0)

p = plot(ts, spacf, size=(1200,800))

end
