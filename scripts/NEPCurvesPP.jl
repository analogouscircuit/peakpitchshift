module NEPCurvesPP

using Plots
using JLD

include("../util/EdgePitchUtils.jl")
include("../util/PeakPick.jl")
include("../util/StimGen.jl")
include("../util/DahlPlot.jl")
include("../core/Sieves.jl")
include("../core/Sigmas.jl")
include("../core/Salience.jl")
include("../core/Envelopes.jl")

using .EdgePitchUtils: edgepitch

# Load HCC Paper Data
hcc_lp = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/hcc_lp_data_clean.jld")
hcc_hp = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/hcc_hp_data_clean.jld")

# Trial Parameters (signals and kernel)
## Load Stimuli
fs = 48000
f0s = [ 200.0 * 2 ^ (k/6) for k ∈ -4:22 ]
dur = 0.1
τ = collect(range(0, stop=dur, step=1/fs))

## Generate SACs
println("Generating stimuli..")
sacs_hp = []
sacs_lp = []
for f0 ∈ f0s
	push!(sacs_hp, EdgePitchUtils.a_hp(τ, f0))
	push!(sacs_lp, EdgePitchUtils.a_lp(τ, f0))
end

## Kernel/Sieve
# σ0 = 4e-4
σ0 = 3e-4

## Max delay parameters
# τ200  = 0.025
# τ2500 = 0.02
τ200  = 0.100
τ2500 = 0.100

# Find Peaks and Period Estimates
println("Finding peaks and making pitch estimates ($(length(sacs_hp)))...")
peaks_hp = []
peaks_lp = []
τ0s_hp = []
τ0s_lp = []
@time for k ∈ 1:length(sacs_hp) 
	print("$k ")
	f0 = f0s[k]
	maxpeaks = Int(floor(Envelopes.hcc_maxdelay(1/f0; τ200=τ200, τ2500=τ2500) * f0 / 2))
	println("$f0 maxpeaks: ", maxpeaks)
	sig(t0) = Sigmas.uniform(σ0=σ0*(155.0/f0), num_bumps=1)
	sieve(t, t0) = Sieves.sievegauss(t, t0, sig(t0); skew=0.075, env=:none)
	peakf(t, t0, x) = Salience.peak_sal(t, t0, x; sievef=sieve)
	sac_hp = sacs_hp[k]
	sac_lp = sacs_lp[k]
	pks_hp = PeakPick.hcpeaks(τ, 1/f0, sac_hp; peakf=peakf)
	pks_lp = PeakPick.hcpeaks(τ, 1/f0, sac_lp; peakf=peakf)
	τ0_hp = EdgePitchUtils.delayfrompeaks(pks_hp; maxpeaks=maxpeaks)
	τ0_lp = EdgePitchUtils.delayfrompeaks(pks_lp; maxpeaks=maxpeaks)
	push!(peaks_hp, pks_hp)
	push!(peaks_lp, pks_lp)
	push!(τ0s_hp, τ0_hp)
	push!(τ0s_lp, τ0_lp)
	if k == length(sacs_hp)
		println()
	end
end

p = plot(xscale=:log)
plot!(p, f0s, (1 ./ τ0s_hp) ./ f0s)
plot!(p, f0s, (1 ./ τ0s_lp) ./ f0s)
scatter!(p, hcc_lp["freqs"], hcc_lp["mt"], yerr=(hcc_lp["yerrlo"], hcc_lp["yerrhi"]))
scatter!(p, hcc_hp["freqs"], hcc_hp["mt"], yerr=(hcc_hp["yerrlo"], hcc_hp["yerrhi"]))

plot!(p, f0s, (edgepitch.(f0s, 0.015; kind=:lp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5)
plot!(p, f0s, (edgepitch.(f0s, 0.015; kind=:hp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5)
plot!(p, f0s, (edgepitch.(f0s, 0.030; kind=:lp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5)
plot!(p, f0s, (edgepitch.(f0s, 0.030; kind=:hp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5)
plot!(p, f0s, (edgepitch.(f0s, 0.060; kind=:lp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5)
plot!(p, f0s, (edgepitch.(f0s, 0.060; kind=:hp)) ./ f0s,
	  linestyle=:dash, linecolor=:black, linewidth=0.5)


end
