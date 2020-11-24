module JASA_NEP_figure_1

using JLD
using Plots
using Measures
using Formatting: format
using Statistics: mean

include("../core/Envelopes.jl")
using .Envelopes: hcc_maxdelay
include("../util/EdgePitchUtils.jl")
using .EdgePitchUtils: edgepitch


################################################################################
# Functions (for statistical analysis)
################################################################################
function lininterp(x, y, x_at)
	@assert length(x) == length(y)
	idx = searchsortedfirst(x, x_at)
	if idx == 1
		return y[idx]
	elseif 1 < idx < length(y)
		return y[idx-1] + (y[idx] - y[idx-1])*(x_at - x[idx-1])/(x[idx]-x[idx-1])
	else
		return y[end]
	end
end

function explainedvariance(μs_obs, μs_model)
	num   = sum( (μs_obs .- μs_model) .^ 2)
	denom = sum( (μs_obs .- mean(μs_obs)) .^ 2)
	100.0 * (1.0 - num/denom)
end

function rmserr(o, p)
	@assert length(o) == length(p)
	(o .- p) .^ 2 |> sum |> sqrt |> x -> x * 100.0 / sqrt(length(o))
end
################################################################################
# Load Data
################################################################################

## Load Data
# Model Data

## Last good ones
# model_lp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_lp_19-11-14-19-36.jld")
# model_hp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_hp_19-11-14-19-19.jld")
# model_lp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_lp_19-11-15-19-4.jld")
# model_hp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_hp_19-11-15-18-47.jld")
# model_lp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_lp_19-11-18-20-5.jld")
# model_hp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_hp_19-11-18-19-18.jld")
# model_lp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_lp_19-11-20-16-1.jld")
# model_hp = load("/home/dahlbom/research/AdaptiveSieves/data/PMPeaksRealNEPCurves_hp_19-11-20-15-13.jld")
## Latest

# model_lp = load("/home/dahlbom/research/AdaptiveSieves/data/latest_lp_good.jld")
# model_hp = load("/home/dahlbom/research/AdaptiveSieves/data/latest_hp_good.jld")

model_lp = load("/home/dahlbom/research/AdaptiveSieves/data/latest_lp.jld")
model_hp = load("/home/dahlbom/research/AdaptiveSieves/data/latest_hp.jld")


## Load HCC data
hcc_lp = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/hcc_lp_data_clean.jld")
hcc_hp = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/hcc_hp_data_clean.jld")


################################################################################
# Plot
################################################################################
xticks=[collect(100.0:100.0:1000.0); 2000.0; 2500.0]
xticklabels = [format("{:.1f}", tick) for tick ∈ xticks]
yticks = 0.80:0.05:1.35
yticklabels = [format("{:.2f}", tick) for tick ∈ yticks]
ylims=(0.78, 1.37)
idcs = [4,6,7,8,9]
xticklabels[idcs] .= ""

lw1 = 2.5
lw2 = 2.0
lw3 = 1.5

p1 = plot(xscale=:log10,
		  xticks=(xticks, xticklabels),
		  ylims=ylims,
		  yticks=(yticks, yticklabels),
		  framestyle=:box,
		  # legend=:none,
		  legendfontsize=12,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  guidefontsize=16,
		  xrotation=45,
		  xlabel="Edge Frequency (Hz)",
		  ylabel="Percent Mistuning")

## Plot HCC Observed Data
plot!(p1, hcc_lp["freqs"], hcc_lp["mt"], yerr=(hcc_lp["yerrlo"], hcc_lp["yerrhi"]),
	  lc=:black,
	  markerstrokecolor=:black,
	  markercolor=:black,
	  label="observed",
	  la=0.75,
	  ls=:dash,
	  lw=lw2)
plot!(p1, hcc_hp["freqs"], hcc_hp["mt"], yerr=(hcc_hp["yerrlo"], hcc_hp["yerrhi"]),
	  lc=:black,
	  markerstrokecolor=:black,
	  markercolor=:black,
	  label="",
	  la=0.75,
	  ls=:dash,
	  lw=lw2)

# Plot HCC Model data
f0s = collect(200*2^(-2/3):0.1:200*2^(11/3))
hcc_fitted_lp = [ edgepitch(fe, hcc_maxdelay(1/fe); kind=:lp)/fe for fe ∈ f0s ]
hcc_fitted_hp = [ edgepitch(fe, hcc_maxdelay(1/fe); kind=:hp)/fe for fe ∈ f0s ]
plot!(p1, f0s, hcc_fitted_lp,
	  label="HCC model",
	  lc=:red,
	  la=0.75,
	  lw=lw1)
plot!(p1, f0s, hcc_fitted_hp,
	  label="",
	  lc=:red,
	  la=0.75,
	  lw=lw1)


# Plot Ideal Curve from fitted data (i.e. non-constant delay)
p2 = plot(xscale=:log10,
		  ylims=ylims,
		  yticks=(yticks, yticklabels),
		  xticks=(xticks, xticklabels),
		  framestyle=:box,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  guidefontsize=16,
		  xrotation=45,
		  xlabel="Edge Frequency (Hz)",
		  ylabel="Percent Mistuning",
		  legendfontsize=12)


## Plot HCC Observed Data
plot!(p2, hcc_lp["freqs"], hcc_lp["mt"], yerr=(hcc_lp["yerrlo"], hcc_lp["yerrhi"]),
	  lc=:black,
	  markerstrokecolor=:black,
	  markercolor=:black,
	  label="observed",
	  la=0.75,
	  ls=:dash,
	  lw=lw2)
plot!(p2, hcc_hp["freqs"], hcc_hp["mt"], yerr=(hcc_hp["yerrlo"], hcc_hp["yerrhi"]),
	  lc=:black,
	  markerstrokecolor=:black,
	  markercolor=:black,
	  label="",
	  la=0.75,
	  ls=:dash,
	  lw=lw2)

# Plot HCC Ideal Curves with constant delays (15, 30, 60 ms)
# plot!(p2, f0s, (edgepitch.(f0s, 0.015; kind=:lp)) ./ f0s,
# 	  linestyle=:dash, linecolor=:black, linewidth=0.5)
# plot!(p2, f0s, (edgepitch.(f0s, 0.015; kind=:hp)) ./ f0s,
# 	  linestyle=:dash, linecolor=:black, linewidth=0.5)
# plot!(p2, f0s, (edgepitch.(f0s, 0.030; kind=:lp)) ./ f0s,
# 	  linestyle=:dash, linecolor=:black, linewidth=0.5)
# plot!(p2, f0s, (edgepitch.(f0s, 0.030; kind=:hp)) ./ f0s,
# 	  linestyle=:dash, linecolor=:black, linewidth=0.5)
# plot!(p2, f0s, (edgepitch.(f0s, 0.060; kind=:lp)) ./ f0s,
# 	  linestyle=:dash, linecolor=:black, linewidth=0.5)
# plot!(p2, f0s, (edgepitch.(f0s, 0.060; kind=:hp)) ./ f0s,
# 	  linestyle=:dash, linecolor=:black, linewidth=0.5)

## Plot first peak data
start_idx = 1
plot!(p2, model_lp["f0s"][start_idx:end], model_lp["mt_μs_peak"][start_idx:end],
	  #yerr=model_lp["mt_stds_peak"][start_idx:end],
	  label="first peak",
	  lc=:blue,
	  la=0.5,
	  ls=:dot,
	  lw=lw2,
	  markercolor=:blue,
	  markerstrokecolor=:blue)
plot!(p2, model_hp["f0s"][start_idx:end], model_hp["mt_μs_peak"][start_idx:end],
	  label="",
	  #yerr=model_hp["mt_stds_peak"][start_idx:end],
	  lc=:blue,
	  la=0.5,
	  ls=:dot,
	  lw=lw2,
	  markercolor=:blue,
	  markerstrokecolor=:blue)

## Plot multi-peak data
# plot!(p2, model_lp["f0s"][start_idx:end], model_lp["mt_μs_hcc"][start_idx:end],
# 	  #yerr=model_lp["mt_stds_hcc"][start_idx:end],
# 	  label="multi-peak",
# 	  lc=:magenta,
# 	  la=0.7,
# 	  ls=:dashdot,
# 	  lw=lw2,
# 	  markercolor=:magenta,
# 	  markerstrokecolor=:magenta)
# plot!(p2, model_hp["f0s"][start_idx:end], model_hp["mt_μs_hcc"][start_idx:end],
# 	  label="",
# 	  #yerr=model_hp["mt_stds_hcc"][start_idx:end],
# 	  lc=:magenta,
# 	  la=0.7,
# 	  ls=:dashdot,
# 	  lw=lw2,
# 	  markercolor=:magenta,
# 	  markerstrokecolor=:magenta)

# Plot model data
plot!(p2, model_lp["f0s"][start_idx:end], model_lp["mt_μs"][start_idx:end],
	  label="multi-peak shift",
	  yerr=model_lp["mt_stds"][start_idx:end],
	  lc=:red,
	  markerstrokecolor=:red,
	  markercolor=:red,
	  lw=lw1)
plot!(p2, model_hp["f0s"][start_idx:end], model_hp["mt_μs"][start_idx:end],
	  label="",
	  yerr=model_hp["mt_stds"][start_idx:end],
	  lc=:red,
	  markerstrokecolor=:red,
	  markercolor=:red,
	  lw=lw1)



p = plot(p2, p1, layout=(2,1), size=(700,1000), margin=5mm)

################################################################################
# Statistical analysis
################################################################################
start_idx = 1
# end_idx = 10
end_idx = length(hcc_lp["freqs"])
f_vals = hcc_lp["freqs"][start_idx:end_idx]

μs_lp_obs = hcc_lp["mt"][start_idx:end_idx]
μs_hp_obs = hcc_hp["mt"][start_idx+1:end_idx+1] 


μs_lp_model = [lininterp(model_lp["f0s"], model_lp["mt_μs"], f) for f ∈ f_vals ]
μs_hp_model = [lininterp(model_hp["f0s"], model_hp["mt_μs"], f) for f ∈ f_vals ]

μs_lp_multipeak  = [lininterp(model_lp["f0s"], model_lp["mt_μs_hcc"], f) for f ∈ f_vals ]
μs_hp_multipeak  = [lininterp(model_hp["f0s"], model_hp["mt_μs_hcc"], f) for f ∈ f_vals ]

μs_lp_peak  = [lininterp(model_lp["f0s"], model_lp["mt_μs_peak"], f) for f ∈ f_vals ]
μs_hp_peak  = [lininterp(model_hp["f0s"], model_hp["mt_μs_peak"], f) for f ∈ f_vals ]

μs_lp_hcc   = [edgepitch(fe, hcc_maxdelay(1/fe); kind=:lp)/fe for fe ∈ f_vals]
μs_hp_hcc   = [edgepitch(fe, hcc_maxdelay(1/fe); kind=:hp)/fe for fe ∈ f_vals]

norm = sqrt(length(μs_lp_obs))

ev_lp_model = explainedvariance((μs_lp_obs), (μs_lp_model))
ev_hp_model = explainedvariance((μs_hp_obs), (μs_hp_model))
ev_model = explainedvariance([μs_lp_obs; μs_hp_obs], [μs_lp_model; μs_hp_model])
rms_model = rmserr([μs_lp_obs; μs_hp_obs], [μs_lp_model; μs_hp_model])
println("EV total, model: $ev_model")
println("RMS total, model: $rms_model\n")

rms_lp_model = rmserr(μs_lp_obs, μs_lp_model)
rms_hp_model = rmserr(μs_hp_obs, μs_hp_model)

ev_lp_hcc   = explainedvariance((μs_lp_obs), (μs_lp_hcc))
ev_hp_hcc   = explainedvariance((μs_hp_obs), (μs_hp_hcc))
ev_hcc = explainedvariance([μs_lp_obs; μs_hp_obs], [μs_lp_hcc; μs_hp_hcc])
rms_hcc = rmserr([μs_lp_obs; μs_hp_obs], [μs_lp_hcc; μs_hp_hcc])
println("EV total, hcc: $ev_hcc")
println("RMS total, hcc: $rms_hcc\n")

rms_lp_hcc = rmserr(μs_lp_obs, μs_lp_hcc)
rms_hp_hcc = rmserr(μs_hp_obs, μs_hp_hcc)

ev_lp_peak  = explainedvariance((μs_lp_obs), (μs_lp_peak))
ev_hp_peak  = explainedvariance((μs_hp_obs), (μs_hp_peak))
ev_peak = explainedvariance([μs_lp_obs; μs_hp_obs], [μs_lp_peak; μs_hp_peak])
rms_peak = rmserr([μs_lp_obs; μs_hp_obs], [μs_lp_peak; μs_hp_peak])
println("EV total, peak: $ev_peak")
println("RMS total, peak: $rms_peak\n")

rms_lp_peak = rmserr(μs_lp_obs, μs_lp_peak) 
rms_hp_peak = rmserr(μs_hp_obs, μs_hp_peak)
 
ev_lp_multipeak  = explainedvariance((μs_lp_obs), (μs_lp_multipeak))
ev_hp_multipeak  = explainedvariance((μs_hp_obs), (μs_hp_multipeak))
ev_multipeak = explainedvariance([μs_lp_obs; μs_hp_obs], [μs_lp_multipeak; μs_hp_multipeak])
rms_multipeak = rmserr([μs_lp_obs; μs_hp_obs], [μs_lp_multipeak; μs_hp_multipeak])
println("EV total, multipeak: $ev_multipeak")
println("RMS total, multipeak: $rms_multipeak\n")

rms_lp_multipeak = rmserr(μs_lp_obs, μs_lp_multipeak) 
rms_hp_multipeak = rmserr(μs_hp_obs, μs_hp_multipeak)

# x as freqs, y as freqs
# f_vals = hcc_lp["freqs"] 
# 
# μs_lp_obs = hcc_lp["mt"] .* f_vals
# μs_hp_obs = hcc_hp["mt"][2:end] .* f_vals
# 
# μs_lp_model = [f*lininterp(model_lp["f0s"], model_lp["mt_μs"], f) for f ∈ f_vals ]
# μs_hp_model = [f*lininterp(model_hp["f0s"], model_hp["mt_μs"], f) for f ∈ f_vals ]
# 
# μs_lp_peak  = [f*lininterp(model_lp["f0s"], model_lp["mt_μs_peak"], f) for f ∈ f_vals ]
# μs_hp_peak  = [f*lininterp(model_hp["f0s"], model_hp["mt_μs_peak"], f) for f ∈ f_vals ]
# 
# μs_lp_hcc   = [fe*edgepitch(fe, hcc_maxdelay(1/fe); kind=:lp)/fe for fe ∈ f_vals]
# μs_hp_hcc   = [fe*edgepitch(fe, hcc_maxdelay(1/fe); kind=:hp)/fe for fe ∈ f_vals]
# 
# ev_lp_model = explainedvariance((μs_lp_obs), (μs_lp_model))
# ev_hp_model = explainedvariance((μs_hp_obs), (μs_hp_model))
# ev_lp_hcc   = explainedvariance((μs_lp_obs), (μs_lp_hcc))
# ev_hp_hcc   = explainedvariance((μs_hp_obs), (μs_hp_hcc))
# ev_lp_peak  = explainedvariance((μs_lp_obs), (μs_lp_peak))
# ev_hp_peak  = explainedvariance((μs_hp_obs), (μs_hp_peak))


println("--------------------")
println("RMS Error, HP, Model: $rms_hp_model")
println("Explained Variance, HP, model: $ev_hp_model")
println()
println("RMS Error, HP, Multi-peak: $rms_hp_multipeak")
println("Explained Variance, HP, Multi-peak: $ev_hp_multipeak")
println()
println("RMS Error, HP, Peak: $rms_hp_peak")
println("Explained Variance, HP, peak: $ev_hp_peak")
println()
println("RMS Error, HP, HCC: $rms_hp_hcc")
println("Explained Variance, HP, hcc: $ev_hp_hcc")
println("--------------------")
println()
println("RMS Error, LP, Model: $rms_lp_model")
println("Explained Variance, LP, model: $ev_lp_model")
println()
println("RMS Error, LP, Multi-peak: $rms_lp_multipeak")
println("Explained Variance, LP, Multi-peak: $ev_lp_multipeak")
println()
println("RMS Error, LP, Peak: $rms_lp_peak")
println("Explained Variance, LP, peaks: $ev_lp_peak")
println()
println("RMS Error, LP, HCC: $rms_lp_hcc")
println("Explained Variance, LP, hcc: $ev_lp_hcc")
println()


end
