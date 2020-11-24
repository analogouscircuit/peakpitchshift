module MTPCurvesSal

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

# Trial Parameters (signals and kernel)
## Load Stimuli
println("Loading data...")
data = load("/home/dahlbom/research/AdaptiveSieves/stimuli/MTPSacs155.jld")
fs = data["fs"]
f0 = data["f0"]
mts = data["mts"]
mtps = data["mtps"]
τ = data["τ"]
sacs_p = data["sacs_p"]

## Kernel/Sieve
σ0 = 4e-4
maxpeaks = Int(floor(Envelopes.hcc_maxdelay(1/f0) * f0))
sig(t0) = Sigmas.uniform(σ0=σ0, num_bumps=maxpeaks)
# sig(t0) = Sigmas.uniform_scaled(t0; σ0=σ0, num_bumps=1)
env(t, t0) = Envelopes.gamma(t, 2.5*t0; peakstart = false)
sieve(t, t0) = Sieves.sievegauss(t, t0, sig(t0); skew=0.075, env=env)
peakf(t, t0, x) = Salience.peak_sal(t, t0, x; sievef=sieve)


# Find Peaks and Period Estimates
println("Finding peaks and making pitch estimates ($(length(mtps)))...")

peaks_p = []
τ0s_p = []
@time for mth ∈ mtps
	print("$mth ")
	global peaks_p
	global τ0s_p
	peaks = []
	τ0s = []
	for sac ∈ sacs_p[mth]
		τ0 = Salience.peak_sal(τ, 1/f0, sac; sievef=sieve)
		push!(τ0s, τ0)
	end
	push!(peaks_p, peaks)
	push!(τ0s_p, τ0s)
end

# Calculate Deviations
## to-do

# Plot Results (Salience Curves, Deviations, and Pitch Curves)
## Salience Curves and Peaks Picked (for sanity check)

## Deviations

## Pitch Curves
println("Plotting results...")
p1 = plot()
for mth ∈ mtps
	global p1
	plot!(p1, mts, 1.0 ./ τ0s_p[mth], label="Partial $(mth)")
end

### For Moore 1985 Comparison
τ0s_avg = sum(τ0s_p) ./ length(τ0s_p)
p2 = plot(mts, 1.0 ./ τ0s_avg)

### For Darwin and Ciocca Comparison
mistuned_f = mts .* 620.0
p3 = plot(mistuned_f, 1.0 ./ τ0s_p[4])
bias = 0.1
dc_curve(Δf; a = 155.0 + bias, k = 0.057, s = 20.0) = a + k*Δf*exp(-(Δf^2)/(2*s^2))
mistuned_f = mts .* 620.0
plot!(p3, mistuned_f, dc_curve.(mistuned_f .- 620.0),
      linestyle=:dash,
      linewidth=0.25,
      color=:gray,
      )

dc_f_points = [550.0, 570.0:10.0:670.0..., 690.0]
dc_f0_points = [154.9, 155.1, 154.75, 154.6, 154.5, 154.65, 155.4, 155.65, 155.75, 155.6, 155.4, 155.25, 155.20]
dc_bar_b = [154.75, 154.9, 154.5, 154.3, 154.25, 154.55, 155.2, 155.5, 155.6, 155.4, 155.25, 155.0, 155.0]
dc_bar_t = [155.1, 155.25, 154.95, 154.8, 154.65, 154.75, 155.0, 155.75, 156.0, 155.75, 155.55, 155.45, 155.45]
err_bottom = [abs(x-y) for (x,y) in (dc_bar_b, dc_f0_points)]
err_top = [abs(x-y) for (x,y) in (dc_bar_t, dc_f0_points)]
#=
err_bottom = [0.1 for k in 1:length(dc_f_points)]
err_top = [0.2 for k in 1:length(dc_f_points)]
=#

scatter!(p3, dc_f_points, dc_f0_points,
         color=:gray,
         yerror = (err_bottom, err_top),
         markershape=:square,
         markersize=8,
         markerstrokewidth=0.5,
         markerstrokecolor=:black,
         markeralpha=0.3,
         markercolor=:black,
         )




p = plot(p1, p2, p3, size=(1200, 400))


end
