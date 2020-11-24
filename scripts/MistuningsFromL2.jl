module MistuningsFromL2

using JLD
using Plots

include("../util/EdgePitchUtils.jl")
include("../util/PeakPick.jl")
using .EdgePitchUtils: delayfrompeaks
using .PeakPick: findpeaknear, peakinterp, hcpeaks

################################################################################
# Script
################################################################################
data = load("/home/dahlbom/research/AdaptiveSieves/stimuli/sac_hc_data_with_fundamental_fft.jld")
sacs = data["sacs"]
fs = data["fs"]
τ = data["delays"]
mistunings = data["mistunings"]

p = plot()
numpeaks = collect(1:2:16)

for maxpeaks ∈ numpeaks
	estimates = Array{Float64, 1}(undef, 0)
	for (n,sac) ∈ enumerate(sacs)
		println("Starting $n of $(length(sacs))")
		peaks = hcpeaks(τ, 1/155.0, sac)
		est = delayfrompeaks(peaks; maxpeaks=maxpeaks)
		push!(estimates, 1/est)
	end

	plot!(p, mistunings .* 620.0, estimates)#ylims=(154.0, 156.5))
end


## Darwin and Ciocca's curve
bias = 0.1
dc_curve(Δf; a = 155.0 + bias, k = 0.057, s = 20.0) = a + k*Δf*exp(-(Δf^2)/(2*s^2))
mistuned_f = mistunings .* 620.0
plot!(p, mistuned_f, dc_curve.(mistuned_f .- 620.0),
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

scatter!(p, dc_f_points, dc_f0_points,
         color=:gray,
         yerror = (err_bottom, err_top),
         markershape=:square,
         markersize=8,
         markerstrokewidth=0.5,
         markerstrokecolor=:black,
         markeralpha=0.3,
         markercolor=:black,
         )

display(p)

end
