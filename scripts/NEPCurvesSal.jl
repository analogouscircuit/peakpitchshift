module NEPCurvesSal

using Plots; pyplot()
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

# Trial Parameters (signals and kernel)
## Load Stimuli
fs = 48000
f0s = [ 200.0 * 2 ^ (k/6) for k ∈ -2:22 ]
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
σ0 = 4e-4

# Find Peaks and Period Estimates
println("Finding peaks and making pitch estimates ($(length(sacs_hp)))...")
τ0s_hp = []
τ0s_lp = []
@time for k ∈ 1:length(sacs_hp) 
	print("$k ")
	f0 = f0s[k]
	maxpeaks = Int(floor(Envelopes.hcc_maxdelay(1/f0) * f0))
	sig(t0) = Sigmas.uniform(σ0=σ0*(155.0/f0), num_bumps=maxpeaks)
	sieve(t, t0) = Sieves.sievegauss(t, t0, sig(t0); skew=0.075, env=:none)
	#peakf(t, t0, x) = Salience.peak_sal(t, t0, x; sievef=sieve)
	sac_hp = sacs_hp[k]
	sac_lp = sacs_lp[k]
	τ0_hp = Salience.peak_sal(τ, 1/f0, sac_hp; sievef=sieve)
	τ0_lp = Salience.peak_sal(τ, 1/f0, sac_lp; sievef=sieve)
	push!(τ0s_hp, τ0_hp)
	push!(τ0s_lp, τ0_lp)
	if k == length(sacs_hp)
		println()
	end
end

p = plot()
plot!(p, f0s, (1 ./ τ0s_hp) ./ f0s)
plot!(p, f0s, (1 ./ τ0s_lp) ./ f0s)

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
