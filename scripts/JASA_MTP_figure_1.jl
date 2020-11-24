module JASA_MTP_figure_1

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
	(o .- p) .^ 2 |> sum |> sqrt |> x -> x / sqrt(length(o))
end
## Load Data
# Model Data

## Last good ones
model_data = load("/home/dahlbom/research/AdaptiveSieves/data/DC_latest.jld")

## Load Darwin-Ciocca Data
dc_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/dc_data_clean.jld")


xaxis = ("Mistuning Factor of 4th Harmonic", font(16))
yaxis = ("Perceived F0 (Hz)", font(16))

p = scatter(dc_data["mistunings"], dc_data["f0s"],
            color=:gray,
            yerror = (dc_data["yerrlo"], dc_data["yerrhi"]),
            markershape=:circle,
            markersize=10,
            markerstrokewidth=0.5,
            markerstrokecolor=:black,
            markeralpha=0.3,
            markercolor=:black,
			legend=:topleft,
			legendfontsize=12,
			label="observed",
			xaxis=xaxis,
			yaxis=yaxis,
			yticks=154.0:0.2:156.1,
			xtickfontsize=14,
			ytickfontsize=14,
			guidefontsize=16,
			size=(700,500),
			framestyle=:box,
			margin=5mm)

# plot!(p, dc_data["mistunings"], dc_data["f0s"], yerr=(dc_data["yerrlo"], dc_data["yerrhi"]),
# 	 lc=:black, la=.3, ls=:dash, legend=:none)

plot!(p, model_data["mt_vals"], model_data["f0s_ests"], lc=:red, label="multi-peak + shift", lw=3.0)
plot!(p, model_data["mt_vals"], model_data["f0_hcc"], lc=:black, label="multi-peak", ls=:dashdot, lw=2.0)
plot!(p, model_data["mt_vals"], model_data["f0_peaks"], lc=:blue, label="first peak", ls=:dot, lw=2.0)


# Calculate Explained Variation
mts_dc = dc_data["mistunings"]
mts_model = model_data["mt_vals"]
μs_dc = dc_data["f0s"]

μs_model = [lininterp(mts_model, model_data["f0s_ests"], mt) for mt ∈ mts_dc]
μs_peak = [lininterp(mts_model, model_data["f0_peaks"], mt) for mt ∈ mts_dc]
μs_hcc = [lininterp(mts_model, model_data["f0_hcc"], mt) for mt ∈ mts_dc]

explained_model = explainedvariance(μs_dc, μs_model) 
rms_model = rmserr(μs_dc, μs_model) 
explained_peak = explainedvariance(μs_dc, μs_peak) 
rms_peak = rmserr(μs_dc, μs_peak) 
explained_hcc = explainedvariance(μs_dc, μs_hcc) 
rms_hcc = rmserr(μs_dc, μs_hcc) 

println("Explained Variation, Model: $explained_model")
println("RMS Error, Model: $rms_model")
println("Explained Variation, Peak: $explained_peak")
println("RMS Error, Peak: $rms_peak")
println("Explained Variation, Multi-peak: $explained_hcc")
println("RMS Error, Multi-peak: $rms_hcc")

end
