module ChannelShiftTest

using AuditoryFilters
using Plots

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/DahlPlot.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/EdgePitchUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Salience.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Sieves.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/Sigmas.jl")
include("/home/dahlbom/research/AdaptiveSieves/scripts/ChannelShiftUtils.jl")
using .ChannelShiftUtils

# -------------------- Parameters --------------------
# Periphery parameters
fs = 50e3
num_chan = 50
f_lo = 100.0
fb = make_erb_filterbank(fs, num_chan, f_lo)
cfs = fb.center_frequencies
cfs_out = exp.( range( log(50.0), log(5000), length=25 ) )
println(cfs_out)

# Stimulus parameters
f0 = 200.0
dur = 0.15
dur_used = 0.05
mth = 4
mtp = 1.060
mistunings = [1.0 for _ ∈ 1:12]
mistunings[mth] = mtp
stim = StimGen.hc(dur, f0, fs; mistunings=mistunings)

# Model Parameters
skew = 0.0
# bump = :gauss
# σ0 = 5.0e-4
bump = :cosbump2
σ0 = 1.5e-3
σ_overall = 2.5e-3 

σ_chan_const = 0.750e-3
σ_sac_const  = 5.0e-3
inn_width = 1.0

const_chan = true

# -------------------- Main Script --------------------
# Filter and process channels
chans = filt(fb, stim)
idx_max = Int(round(dur_used*fs))
println("max idx: $idx_max")
chans = chans[1:idx_max,:]
ts = collect(0:idx_max-1) ./ fs

# Half-wave rectify
chans = map( x -> x > 0.0 ? x : 0.0, chans)
chans_hwr = copy(chans)

# Autocorrelate
chans = ac_channels(chans)
chans_ac = copy(chans)
peaks_cf_ac = Float64[]
for n ∈ 1:size(chans_ac)[2]
	peaks = PitchModelPeaks.findpeaks(ts, chans_ac[:,n], 1/f0, fs, 1)
	push!(peaks_cf_ac, peaks[1])
end

# Kill first bump
# for n ∈ 1:size(chans)[2]
# 	cf = cfs[n]
# 	t_max = 1.0/(2.0*cf)
# 	idx_max = Int(round(t_max*fs))
# 	chans[1:idx_max,n] .= 0.0
# end

# Local summation
weightf(cf_out, cf_in) = w_gauss_full(cf_out, cf_in; σ=inn_width)
chans = sum_channels_locally(chans, cfs, cfs, weightf)
chans_summed = copy(chans)

# Smooth locally summed channels
#sigmaf = _ -> σ0 * 155.0/f0
# sigmaf = _ -> 0.0001
# sigmaf(cf) = 4.0 * σ0 * 155.0 / cf 	# determine smoothing kernel width
# sigmaf(cf) = 4.0 * 3.0e-4 * 155.0 / cf
if const_chan
	sigmaf(cf) = σ_chan_const 
else
	sigmaf(cf) = σ0 * 200.0 / cf
end
chans = smooth_channels(chans, cfs, sigmaf, fs)
chans ./= maximum(chans)/5.0
chans_smoothed = copy(chans)

# Generate final SACF
sac = sum(chans, dims=2)[:,1]
#sig(t0) = Sigmas.uniform(σ0=σ_overall, num_bumps=1)
sig(t0) = Sigmas.uniform(σ0=σ_sac_const, num_bumps=1)
sieve(t, t0) = Sieves.sievefunc(t, t0, sig(t0); bump=bump, skew=skew, env=:none)
sac_smoothed = Salience.salience(ts, ts, sac, sieve)
sac = sac_smoothed


## Plot results
xlims=(0.000, 0.025)
fill_attr = (0, 1.0, :white)
p_chans = DahlPlot.stagger_plot_mat(ts, chans;
									offset=0.05, rev=true, xlims=(0.0, 0.05))
p_hwr = DahlPlot.stagger_plot_mat(ts, chans_hwr;
								  offset=0.05, rev=true, xlims=(0.0, 0.05))


p_ac = DahlPlot.stagger_plot_mat(ts, chans_ac;
								 offset=0.5, rev=true, xlims=xlims)
plot!(p_ac, [1/f0, 1/f0], [0, 25.0], lc=:red, ls=:dash, title="ac")
p_ac_sac = plot(ts, sum(chans_ac, dims=2)[:,1], color=:black, xlims=xlims)
plot!(p_ac_sac, [1/f0, 1/f0], [0, 30], color=:red, ls=:dash)
l = @layout [a{0.8h}; b]
p_ac = plot(p_ac, p_ac_sac, layout=l)



p_summed = DahlPlot.stagger_plot_mat(ts, chans_summed;
									 offset=0.5, rev=true, xlims=xlims)
plot!(p_summed, [1/f0, 1/f0], [0, 25.0], lc=:red, ls=:dash)
p_summed_sac = plot(ts, sum(chans_summed, dims=2)[:,1], color=:black, xlims=xlims)
plot!(p_summed_sac, [1/f0, 1/f0], [0, 300], color=:red, ls=:dash)
l = @layout [a{0.8h}; b]
p_summed = plot(p_summed, p_summed_sac, layout=l)


p_smoothed = DahlPlot.stagger_plot_mat(ts, chans_smoothed;
									   offset=0.5, rev=true, xlims=xlims)
plot!(p_smoothed, [1/f0, 1/f0], [0, 25.0], lc=:red, ls=:dash)
#p_smoothed_sac = plot(ts, sum(chans_smoothed, dims=2)[:,1], color=:black, xlims=xlims)
p_smoothed_sac = plot(ts, sac, color=:black, xlims=xlims)
plot!(p_smoothed_sac, [1/f0, 1/f0], [0, 150], color=:red, ls=:dash)
l = @layout [a{0.8h}; b]
p_smoothed = plot(p_smoothed, p_smoothed_sac, layout=l)


p = plot(p_ac, p_summed, p_smoothed, layout=(1,3), size=(1500,800))

p_sac = plot(ts, sac)
maxpeaks = 3
peaks = PitchModelPeaks.findpeaks(ts, sac, 1/f0, fs, maxpeaks)
τ0 = EdgePitchUtils.delayfrompeaks(peaks; maxpeaks=maxpeaks)
println("Estimate: ", 1/τ0)

end
