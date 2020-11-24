module IdealDeviations

using Plots
include("../util/PeakPick.jl")
include("../core/Envelopes.jl")

"""
	hc(t, f0, n)

Harmonic complex with fundamental of f0 and n harmonics.
"""
function hc(t, f0, n)
	0.5 .* ( sin.( (n + 0.5) * 2π * f0 .* t ) ./ sin.( 0.5 * 2π * f0 .* t )  .- 1.0 )	
	#=
	out = zeros(length(t))
	for k ∈ 1:n
		out .+= cos.( 2π * f0 * k .* t)
	end
	return out
	=#
end


"""
	dhc(t, f0, n)

Derivative of a harmonic complex with a fundamental of f0 and n harmonics.
"""
function dhc(t, f0, n)
	a1(x) = π * (n + 0.5) * f0 * cos( (n + 0.5) * 2 * π * f0 * x ) / sin(π * f0 * x)
	a2(x) = - π * f0 * sin( (n + 0.5) * 2 * π * f0 * x ) * cos( π * f0 * x ) / 
			  ( 2 * sin( π * f0 * x) ^ 2)
	d = a1.(t) + a2.(t)
	map( x -> isnan(x) || x == Inf ? 0.0 : x, d)
	#=
	out = zeros(length(t))
	for k ∈ 1:n
		out .+= k .* sin.(2π * f0 * k .* t)
	end
	return -2π .* out
	=#
end


"""
	hcmt(t, f0, n, m, p)

Harmonic complex with a fundamental of f0 and n harmonics.  The m-th harmonic is mistuned by a factor of p.
"""
function hcmt(t, f0, n, m, p)
	h(t, f0, n) .+ cos.(2π * m * p * f0 .* t) .- cos.( 2π * m * f0 .* t)
end

"""
	dhcmt(t, f0, n, m, p)

Derivative of a harmonic complex with a mistuned partial.  Fundamental of f0, n harmonics, m-th harmonic detuned by factor of p.
"""
function dhcmt(t, f0, n, m, p)
	#=
	a1(x) = π * (n + 0.5) * f0 * cos( (n + 0.5) * 2 * π * f0 * x ) / sin(π * f0 * x)
	a2(x) = - π * f0 * sin( (n + 0.5) * 2 * π * f0 * x ) * cos( π * f0 * x ) / 
			  ( 2 * sin( π * f0 * x) ^ 2)
	a3(x) = 2 * π * m * f0 * ( sin( 2π * m * f0 * x ) - p * sin(2π * m * p * f0 * x) )
	return a1.(t) + a2.(t) + a3.(t)
	=#
	out = zeros(length(t))
	for k ∈ 1:n
		out .+= k .* sin.( 2π * k * f0 .* t)
	end
	out .*= -2π * f0
	out .+= 2π * m * f0 .* sin.(2π * m * f0 .* t)
	out .-= 2π * m * p * f0 .* sin.(2π * m * p *f0 .* t)
	if length(out) == 1
		return out[1]
	end
	return out
end


"""
"""
function ddhcmt(t, f0, n, m, p)
	out = zeros(length(t))
	for k ∈ 1:n
		out .+= k^2 .* cos.( 2π * k * f0 .* t)
	end
	out .*= - (2π * f0) ^ 2
	out .+= ((2π * m * f0) ^ 2) .* cos.(2π * m * f0 .* t)
	out .-= ((2π * m * p * f0) ^ 2) .* cos.(2π * m * p * f0 .* t)
	if length(out) == 1
		return out[1]
	end
	return out
end


"""
	newton(f, df, x0; tol=1e-10, maxits = 100)

Newston-Raphson root finding method.
"""
function newton(f, df, x0; tol=1e-9, maxits = 100, verbose=true)
	count = 0
	x_new = 0.0
	while count < maxits 
		x_new = x0 - f(x0) / df(x0)
		if abs(x_new - x0) < tol
			break
		end
		count += 1
		x0 = x_new
	end
	if verbose
		if count == maxits
			println("Didn't converge!")
		end
	end
	return x_new
end


function findpeaks(t_max, f0, f, df)
	num_peaks = floor(t_max * f0)
	t0s = [k/f0 for k ∈ 1:num_peaks]
	peaks = []
	for t0 in t0s
		push!(peaks, newton(f, df, t0))
	end
	peaks
end


# Main Script
fs = 100000
dur = 0.1
t = collect(range(0, stop=dur, step=1/fs)) 
f0s = [200 * 2 ^ k for k in 1:5]
n = 12
m = 4
p = 0.965 
ps = [ 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06]

plots_dev = []
plots_con = []
plots_dev_rel = []
mags = []
#for f0 in f0s
f0 = 400.0
for p in ps
	println("$f0")
	n_peaks = floor(dur*f0)
	peaks = [k/f0 for k in 1:n_peaks]
	w = Envelopes.gamma(peaks, 1/f0; peakstart=false)
	f(x)  = dhcmt(x, f0, n, m, p)
	df(x) = ddhcmt(x, f0, n, m, p)
	peaks_hcmt = findpeaks(t[end], f0, f, df) 
	devs = peaks .- peaks_hcmt
	devs .* f0
	p_dev = plot(devs, line=:stem, marker=:circle, legend=:none, title="$p")
	devs_rel = devs .* [1/k for k ∈ 1:length(devs)]
	devs_rel = (1.0 .- ( 1.0 ./ ((1/f0) .+ devs_rel)) ./ (f0)) 
	#devs_rel = 1.0 .- ((1.0/f0) .+ devs_rel) ./ (1.0 / f0)
	p_dev_rel = plot(devs_rel, line=:stem, marker=:circle, legend=:none, title="$p")#, ylims=(-2.5e-6, 2.5e-6))
	p_con = plot(devs .* w, ylims=(-2.5e-6, 2.5e-6), line=:stem, marker=:circle, legend=:none, title="$p")
	push!(plots_dev, p_dev)
	push!(plots_con, p_con)
	push!(plots_dev_rel, p_dev_rel)
	push!(mags, sum(devs .* w))
end

p_devs = plot(plots_dev...)
p_rel  = plot(plots_dev_rel...)
p_cons = plot(plots_con...)
p_tots = plot(ps, mags)





end
