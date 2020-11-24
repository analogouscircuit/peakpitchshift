module JASAMooreDataFigPltMR

using Plots 
using JLD
using Statistics: mean, std

moore_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/moore_complete_data.jld")
model_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/model_moore_complete.jld")

moore_data = moore_data["data_all"]
model_data = model_data["data_all"]

subjects = [:BG, :BM, :MS]
f0s = [100, 200, 400]
partials = 1:5  # disregard sixth
mts = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0]
max_model_p = 6
max_moore_p = 5

################################################################################
# By partial 
################################################################################
# average over frequencies -- moore data
subject_partials = [[] for _ ∈ 1:3]
for (k, subject) ∈ enumerate(keys(moore_data))
	global subjects_partial
	partial_pts = []
	for partial ∈ 1:max_moore_p
		for f0 ∈ f0s
			push!(partial_pts, moore_data[subject][f0][partial])
		end
		push!(subject_partials[k], mean(partial_pts))
	end
end

# average over frequencies -- model data
model_partials = []
for partial ∈ 1:max_model_p
	global model_partials
	partial_pts = []
	for f0 ∈ f0s
		push!(partial_pts, model_data[f0][partial])
	end
	push!(model_partials, mean(partial_pts))
end

# sort by median-range
moore_partials = []
for partial ∈ 1:max_moore_p
	global moore_partials
	p_med = []
	p_lo = []
	p_hi = []
	for mtn ∈ 1:6
		vals = []
		for subject ∈ 1:3
			push!(vals, subject_partials[subject][partial][mtn])
		end
		vals = sort(vals)
		(lo, med, hi) = (vals[1], vals[2], vals[3])
		push!(p_lo, lo)
		push!(p_med, med)
		push!(p_hi, hi)
	end
	push!(moore_partials, [p_lo, p_med, p_hi])
end

# plot compared to model
plts_p = []
for partial ∈ 1:max_moore_p
	p = scatter(mts, moore_partials[partial][2],
		 yerr=(moore_partials[partial][1], moore_partials[partial][3]),
		 ylims=(0.0, 1.75),
		 label="Moore et al.")
	plot!(p, mts, model_partials[partial],
		  label="Model")
	push!(plts_p, p)
end
p_partials = plot(plts_p...)


################################################################################
# By Frequency 
################################################################################


end
