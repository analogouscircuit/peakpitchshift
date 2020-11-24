module JASAMooreDataFigPlot

using Plots 
using Measures
using JLD
using Statistics: mean, std

moore_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/moore_complete_data.jld")
model_data = load("/home/dahlbom/research/AdaptiveSieves/data/paper_data/model_moore_complete.jld")

moore_data = moore_data["data_all"]
peak_data = model_data["data_all_peak"]
model_data = model_data["data_all"]

subjects = [:BG, :BM, :MS]
f0s = [100, 200, 400]
partials = 1:5  # disregard sixth
mts = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0]


################################################################################
# By f0
################################################################################
max_model_p = 6
ylims_f=(-0.15, 1.3)
ylims_p=(-0.1, 1.70)
ylims_all=(-0.0, 0.8)

# -------------------- 100 Hz --------------------
model_100 = zeros(length(mts))
model_100_peaks = zeros(length(mts))
count = 0
for k ∈ 1:max_model_p
	global model_100 += model_data[100][k]
	global model_100_peaks += peak_data[100][k]
	global count += 1
end
model_100 ./= count
model_100_peaks ./= count

subject_100 = zeros(length(mts))
s100_points = []
count = 0
for subject ∈ keys(moore_data)
	for k ∈ 1:5
		# global subject_100 += moore_data[subject][100][k]
		global count += 1
		global s100_points
		push!(s100_points, moore_data[subject][100][k])
	end
end
subject_100_μ = mean(s100_points) 
subject_100_std = std(s100_points)

p100 = scatter(mts, subject_100_μ, yerr=subject_100_std, label="moore", ylims=ylims_f, legend=:none,
			   color=:gray,
			   # markershape=:circle,
			   markersize=6,
			   markerstrokewidth=0.5,
			   markerstrokecolor=:black,
			   markeralpha=0.3,
			   markercolor=:black,
			   title="100 Hz",
			   xlabel="Partial Mistuning (%)",
			   ylabel="Pitch Mistuning (%)",
			   framestyle=:box)
plot!(p100, mts, model_100, label="model", lc=:red, lw=2.5)
plot!(p100, mts, model_100_peaks, label="first peak", lc=:blue, ls=:dot, lw=2.0)


# -------------------- 200 Hz --------------------
model_200 = zeros(length(mts))
model_200_peaks = zeros(length(mts))
count = 0
for k ∈ 1:max_model_p
	global model_200 += model_data[200][k]
	global model_200_peaks += peak_data[200][k]
	global count += 1
end
model_200 ./= count
model_200_peaks ./= count

subject_200 = zeros(length(mts))
s200_points = []
count = 0
for subject ∈ keys(moore_data)
	for k ∈ 1:5
		# global subject_200 += moore_data[subject][200][k]
		global count += 1
		global s200_points
		push!(s200_points, moore_data[subject][200][k])
	end
end
subject_200_μ = mean(s200_points) 
subject_200_std = std(s200_points)

p200 = scatter(mts, subject_200_μ, yerr=subject_200_std, label="moore", ylims=ylims_f, legend=:none,
			   color=:gray,
			   markershape=:circle,
			   markersize=6,
			   markerstrokewidth=0.5,
			   markerstrokecolor=:black,
			   markeralpha=0.3,
			   markercolor=:black,
			   title="200 Hz",
			   xlabel="Partial Mistuning (%)",
			   framestyle=:box)
plot!(p200, mts, model_200, label="model", lc=:red, lw=2.5)
plot!(p200, mts, model_200_peaks, label="first peak", lc=:blue, ls=:dot, lw=2.0)


# -------------------- 400 Hz --------------------
model_400 = zeros(length(mts))
model_400_peaks = zeros(length(mts))
count = 0
for k ∈ 1:max_model_p
	global model_400 += model_data[400][k]
	global model_400_peaks += peak_data[400][k]
	global count += 1
end
model_400 ./= count
model_400_peaks ./= count

subject_400 = zeros(length(mts))
s400_points = []
count = 0
for subject ∈ keys(moore_data)
	for k ∈ 1:5
		# global subject_400 += moore_data[subject][400][k]
		global count += 1
		global s400_points
		push!(s400_points, moore_data[subject][400][k])
	end
end
subject_400_μ = mean(s400_points) 
subject_400_std = std(s400_points)

p400 = scatter(mts, subject_400_μ, yerr=subject_400_std, label="moore", ylims=ylims_f, legend=:none,
			   color=:gray,
			   markershape=:circle,
			   markersize=6,
			   markerstrokewidth=0.5,
			   markerstrokecolor=:black,
			   markeralpha=0.3,
			   markercolor=:black,
			   title="400 Hz",
			   xlabel="Partial Mistuning (%)",
			   framestyle=:box)
plot!(p400, mts, model_400, label="model", lc=:red, lw=2.5)
plot!(p400, mts, model_400_peaks, label="first peak", ls=:dot, lc=:blue, lw=2.0)

# -------------------- Final plot by f0 --------------------
p1 = plot(p100, p200, p400, layout=(1,3))


################################################################################
# All
################################################################################

subject_all_μ = mean([subject_400_μ, subject_200_μ, subject_100_μ])
subject_all_std = std([subject_400_μ, subject_200_μ, subject_100_μ])
model_all = mean([model_100, model_200, model_400])

p2 = scatter(mts, subject_all_μ, label="moore", yerr=subject_all_std, ylims=ylims_all, legend=:none)
plot!(p2, mts, model_all, label="model")


################################################################################
# By partial 
################################################################################

# -------------------- Partial 1 --------------------
subject_p1_points = []
for subject ∈ keys(moore_data)
	for f0 ∈ keys(moore_data[subject])
		global subject_p1_points
		push!(subject_p1_points, moore_data[subject][f0][1])
	end
end
subject_p1_μ = mean(subject_p1_points)
subject_p1_std = std(subject_p1_points)

model_p1 = []
model_p1_peaks = []
for f0 ∈ keys(model_data)
	global model_p1
	global model_p1_peaks
	push!(model_p1, model_data[f0][1])
	push!(model_p1_peaks, peak_data[f0][1])
end
model_p1_μ = mean(model_p1)
model_p1_peaks_μ = mean(model_p1_peaks)

pp1 = scatter(mts, subject_p1_μ, yerr=subject_p1_std, ylims=ylims_p, legend=:none,
		   color=:gray,
		   markershape=:circle,
		   markersize=6,
		   markerstrokewidth=0.5,
		   markerstrokecolor=:black,
		   markeralpha=0.3,
		   markercolor=:black,
		   title="1st Harmonic",
		   ylabel="Pitch MT (%)",
		   xlabel="Partial Mistuning (%)",
		   framestyle=:box)
plot!(pp1, mts, model_p1_μ, lc=:red, lw=2.5)
plot!(pp1, mts, model_p1_peaks_μ, lc=:blue, ls=:dot, lw=2.0)

# -------------------- Partial 2 --------------------
subject_p2_points = []
for subject ∈ keys(moore_data)
	for f0 ∈ keys(moore_data[subject])
		global subject_p2_points
		push!(subject_p2_points, moore_data[subject][f0][2])
	end
end
subject_p2_μ = mean(subject_p2_points)
subject_p2_std = std(subject_p2_points)

model_p2 = []
model_p2_peaks = []
for f0 ∈ keys(model_data)
	global model_p2
	global model_p2_peaks
	push!(model_p2, model_data[f0][2])
	push!(model_p2_peaks, peak_data[f0][2])
end
model_p2_μ = mean(model_p2)
model_p2_peaks_μ = mean(model_p2_peaks)

pp2 = scatter(mts, subject_p2_μ, yerr=subject_p2_std, ylims=ylims_p, legend=:none,
		   color=:gray,
		   markershape=:circle,
		   markersize=6,
		   markerstrokewidth=0.5,
		   markerstrokecolor=:black,
		   markeralpha=0.3,
		   markercolor=:black,
		   title="2nd Harmonic",
		   xlabel="Partial Mistuning (%)",
		   framestyle=:box)
plot!(pp2, mts, model_p2_μ, lc=:red, lw=2.5)
plot!(pp2, mts, model_p2_peaks_μ, lc=:blue, ls=:dot, lw=2.0)

# -------------------- Partial 3 --------------------
subject_p3_points = []
for subject ∈ keys(moore_data)
	for f0 ∈ keys(moore_data[subject])
		global subject_p3_points
		push!(subject_p3_points, moore_data[subject][f0][3])
	end
end
subject_p3_μ = mean(subject_p3_points)
subject_p3_std = std(subject_p3_points)

model_p3 = []
model_p3_peaks = []
for f0 ∈ keys(model_data)
	global model_p3
	global model_p3_peaks
	push!(model_p3, model_data[f0][3])
	push!(model_p3_peaks, peak_data[f0][3])
end
model_p3_μ = mean(model_p3)
model_p3_peaks_μ = mean(model_p3_peaks)

pp3 = scatter(mts, subject_p3_μ, yerr=subject_p3_std, ylims=ylims_p, legend=:none,
		   color=:gray,
		   markershape=:circle,
		   markersize=6,
		   markerstrokewidth=0.5,
		   markerstrokecolor=:black,
		   markeralpha=0.3,
		   markercolor=:black,
		   title="3rd Harmonic",
		   xlabel="Partial Mistuning (%)",
		   framestyle=:box)
plot!(pp3, mts, model_p3_μ, lc=:red, lw=2.5)
plot!(pp3, mts, model_p3_peaks_μ, lc=:blue, ls=:dot, lw=2.0)


# -------------------- Partial 4 --------------------
subject_p4_points = []
for subject ∈ keys(moore_data)
	for f0 ∈ keys(moore_data[subject])
		global subject_p4_points
		push!(subject_p4_points, moore_data[subject][f0][4])
	end
end
subject_p4_μ = mean(subject_p4_points)
subject_p4_std = std(subject_p4_points)

model_p4 = []
model_p4_peaks = []
for f0 ∈ keys(model_data)
	global model_p4
	global model_p4_peaks
	push!(model_p4, model_data[f0][4])
	push!(model_p4_peaks, peak_data[f0][4])
end
model_p4_μ = mean(model_p4)
model_p4_peaks_μ = mean(model_p4_peaks)

pp4 = scatter(mts, subject_p4_μ, yerr=subject_p4_std, ylims=ylims_p, legend=:none,
		   color=:gray,
		   markershape=:circle,
		   markersize=6,
		   markerstrokewidth=0.5,
		   markerstrokecolor=:black,
		   markeralpha=0.3,
		   markercolor=:black,
		   title="4th Harmonic",
		   ylabel="Pitch MT (%)",
		   xlabel="Partial Mistuning (%)",
		   framestyle=:box)
plot!(pp4, mts, model_p4_μ, lc=:red, lw=2.5)
plot!(pp4, mts, model_p4_peaks_μ, lc=:blue, ls=:dot, lw=2.0)

# -------------------- Partial 5 --------------------
subject_p5_points = []
for subject ∈ keys(moore_data)
	for f0 ∈ keys(moore_data[subject])
		global subject_p5_points
		push!(subject_p5_points, moore_data[subject][f0][5])
	end
end
subject_p5_μ = mean(subject_p5_points)
subject_p5_std = std(subject_p5_points)

model_p5 = []
model_p5_peaks = []
for f0 ∈ keys(model_data)
	global model_p5
	global model_p5_peaks
	push!(model_p5, model_data[f0][5])
	push!(model_p5_peaks, peak_data[f0][5])
end
model_p5_μ = mean(model_p5)
model_p5_peaks_μ = mean(model_p5_peaks)

pp5 = scatter(mts, subject_p5_μ, yerr=subject_p5_std,
		   color=:gray,
		   markershape=:circle,
		   markersize=6,
		   markerstrokewidth=0.5,
		   markerstrokecolor=:black,
		   markeralpha=0.3,
		   markercolor=:black,
		   title="5th Harmonic",
		   xlabel="Partial Mistuning (%)",
		   framestyle=:box)
plot!(pp5, mts, model_p5_μ, lw=2.5, lc=:red)
plot!(pp5, mts, model_p5_peaks_μ, lc=:blue, ls=:dot, ylims=ylims_p, legend=:none, lw=2.0)

pp6 = scatter(mts, subject_p5_μ, yerr=subject_p5_std,
			  legend=:bottom,
			  legendfontsize=10,
			  label="observed",
		   	  color=:gray,
		   	  markershape=:circle,
		   	  markersize=6,
		   	  markerstrokewidth=0.5,
		   	  markerstrokecolor=:black,
		   	  markeralpha=0.3,
		   	  markercolor=:black,
			  xlims=(10,20),
			  grid=:none,
		   	  framestyle=:none)
plot!(pp6, mts, model_p5_μ, lw=2.5, lc=:red, label="multi-peak+shift")
plot!(pp6, mts, model_p5_peaks_μ, lc=:blue, ls=:dot, ylims=ylims_p, lw=2.0,
	  label="first peak")


p3 = plot(pp1, pp2, pp3, pp4, pp5, pp6, layout=(2,3))

l = @layout [ a{0.333h} ; b ] 

p = plot(p1, p3, layout=l, size=(800,590), margin=4mm)

end
