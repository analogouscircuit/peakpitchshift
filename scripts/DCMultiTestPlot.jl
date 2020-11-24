module DCMultiTestPlot

using JLD
using Plots
using Formatting: format

include("../core/Envelopes.jl")
using .Envelopes: hcc_maxdelay
include("../util/EdgePitchUtils.jl")
using .EdgePitchUtils: edgepitch



## Load Data
# Model Data

## Last good ones
model_data = load("/home/dahlbom/research/AdaptiveSieves/data/DCMulti_latest.jld")
mt_vals = model_data["mt_vals"]
println("mt_vals: $mt_vals")
f0_ests = model_data["f0_ests"]
println("Type of f0_ests: ", typeof(f0_ests))
mt_nums = sort(collect(keys(f0_ests)))
println(mt_nums)


## Load Darwin-Ciocca Data
dc_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/dc_data_clean.jld")


# xaxis = ("Mistuning Factor of 4h Harmonic", font(16))
# yaxis = ("Perceived F0 (Hz)", font(16))
# 
# p1 = scatter(dc_data["mistunings"], dc_data["f0s"],
#             color=:gray,
#             yerror = (dc_data["yerrlo"], dc_data["yerrhi"]),
#             markershape=:square,
#             markersize=8,
#             markerstrokewidth=0.5,
#             markerstrokecolor=:black,
#             markeralpha=0.3,
#             markercolor=:black,
# 			legend=:none,
# 			xaxis=xaxis,
# 			yaxis=yaxis,
# 			yticks=154.0:0.2:156.1,
# 			xtickfontsize=14,
# 			ytickfontsize=14,
# 			guidefontsize=16,
# 			size=(700,500),
# 			framestyle=:box)

# plot!(p1, dc_data["mistunings"], dc_data["f0s"], yerr=(dc_data["yerrlo"], dc_data["yerrhi"]),
# 	 lc=:black, la=.3, ls=:dash, legend=:none)
# 
# plot!(p1, mt_vals, f0_ests[1], lc=:red, legend=:none, lw=2.0)

plots_other = []
for k ∈ mt_nums 
	xaxis = ("Partial $(mt_nums[k])", font(12))
	yaxis = ("Perceived F0 (Hz)", font(12))
	pt = scatter(dc_data["mistunings"], dc_data["f0s"],
				 color=:gray,
				 yerror = (dc_data["yerrlo"], dc_data["yerrhi"]),
				 markershape=:square,
				 markersize=8,
				 markerstrokewidth=0.5,
				 markerstrokecolor=:black,
				 markeralpha=0.3,
				 markercolor=:black,
				 legend=:none,
				 xaxis=xaxis,
				 yaxis=yaxis,
				 yticks=154.0:0.2:156.1,
				 xtickfontsize=14,
				 ytickfontsize=14,
				 guidefontsize=16,
				 size=(700,500),
				 title="Partial Num: $k",
				 framestyle=:box)
	plot!(pt, dc_data["mistunings"], dc_data["f0s"], yerr=(dc_data["yerrlo"], dc_data["yerrhi"]),
		  lc=:black, la=.3, ls=:dash, legend=:none)
	plot!(pt, mt_vals, f0_ests[k], lc=:red, legend=:none, lw=2.0)
	push!(plots_other, pt)
end

# p2 = plot(plots_other..., layout=(2,2))
# 
# p = plot(p1, p2, layout=(1,2), size=(1500,800))
p = plot(plots_other...)

end
