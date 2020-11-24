module PeakDeviations

using Plots

include("../util/PeakPick.jl")
using .PeakPick: peakinterp, hcpeaks

function devplot(τ, τ0, sac; tol=0.15, ylims=(-0.5, 0.5))
	peaks = hcpeaks(τ, τ0, sac; tol=tol)
	p1 = plot(τ, sac, linealpha=0.65, legend=:none)
	diffs = []
	for (k, peak) ∈ enumerate(peaks)
		plot!(p1, [k*τ0, k*τ0], [0.0, 1.0],
			  color=:black, linewidth=0.75, legend=:none)
		plot!(p1, [peak, peak], [0.0, 1.0],
			  color=:red, linestyle=:dash, linewidth=0.75, legend=:none)
		push!(diffs, ((1.0/peak) - (1.0/(k*τ0)))/((1.0/(k*τ0))) )
	end
	p2 = plot(diffs * 100.0, line=:stem, marker=:circle, ylims=ylims, legend=:none)
	# l = @layout [a{0.8h}; b]
	# p = plot(p1, p2, layout=l)
	return p1, p2
end


end
