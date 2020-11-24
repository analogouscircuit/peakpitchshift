module DahlPlot

using Plots

export stagger_plot, stagger_plot_mat, plotpeaks

function stagger_plot(x, y; offset=0.5, xlims=(0,0.04))
      fill_attr = (0, 1.0, :white)
      p = plot(x, y[end],
                  xlims=xlims,
                  legend=:none,
                  color=:black,
                  linewidth=0.5,
                  fill = fill_attr)
      for k ∈ length(y)-1:-1:1
            plot!(x, y[k] .+ (k-1)*offset,
                  color=:black,
                  linewidth=0.65,
                  fill = fill_attr)
      end
      return p
end

function stagger_plot_mat(x, yMat; offset=0.5, rev=false, xlims=(0.0, 0.1))
	if rev
		yMat = reverse(yMat, dims=2)
	end
      fill_attr = (0, 1.0, :white)
      p = plot(x, yMat[:,end],
			   	  xlims=xlims,
                  legend=:none,
                  color=:black,
                  linewidth=0.5,
                  fill = fill_attr)
      for k ∈ size(yMat)[2]-1:-1:1
            plot!(x, yMat[:,k] .+ (k-1)*offset,
                  color=:black,
                  linewidth=0.65,
                  fill = fill_attr)
      end
      return p
end

function plotpeaks(plt, peaks, height; color=:black)
    for peak ∈ peaks
        plot!(plt, [peak, peak], [0, height]; color=color, linewidth=0.64, linestyle=:dash, legend=:none)
    end
    return plt
end

end
