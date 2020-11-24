module EdgePitchCurves

using JLD
using Plots

include("../core/Adapt.jl")
include("../core/Envelopes.jl")
include("../core/Sigmas.jl")
include("../core/Sieves.jl")
include("../util/EdgePitchUtils.jl")

using .EdgePitchUtils: edgepitch, a_lp, a_hp

################################################################################
#
################################################################################


################################################################################
# Main Script
################################################################################
data = load("/home/dahlbom/research/AdaptiveSieves/stimuli/sac_edge_data_avg.jld")
fs = data["fs"]
f0s = data["f0s"]
sacs_lp = data["sacs_lp"]
sacs_hp = data["sacs_hp"]
τ = data["τ"]

function μ_est(f)
	# for σ0 = 2.5e-4, gamma, 20 bumps, uniform_scaled
	idcs = [1, 6, 11, 16, 21, 25]
	f_vals = f0s[idcs]
	μ_vals = [5e-8, 2e-8, 9e-9, 4e-9, 2e-9, 9e-10]
	if f < f_vals[1]
		μ = μ_vals[1]
		return μ
	elseif f >= f_vals[end] 
		μ = μ_vals[end]
		return μ
	end
	for k in 1:(length(f_vals)-1)
		if (f_vals[k] <= f) && (f < f_vals[k+1])
			c = (f - f_vals[k]) / (f_vals[k+1] - f_vals[k])
			μ = c * μ_vals[k] + (1.0 - c) * μ_vals[k+1]
			return μ
		end
	end
end

σ0 = 2.5e-4 #width at τ0 = 1/155
sig(τ0) = Sigmas.uniform_scaled(τ0; σ0=σ0, num_bumps=20)
#sig(τ0) = Sigmas.nonuniform_scaled(τ0; σ0=σ0, num_bumps=5)
#sig(τ0) = Sigmas.uniform(; σ0=σ0, num_bumps=5)
# env = Envelopes.hcc
env = Envelopes.gamma
#sig(τ0) = Sigmas.nonuniform_envinv(σ0, env; σ0=0.0, num_bumps=5)
sieve(τ, τ0) = Sieves.sievegauss(τ, τ0, sig(τ0); env=env)

p1 = plot(τ, sieve(τ, 1/200))
plot!(τ, sieve(τ, 1/400))
plot!(τ, sieve(τ, 1/800))
plot!(τ, sieve(τ, 1/1600))

ests_lp = []
ests_hp = []
for k in 1:length(sacs_hp)
	println("$(f0s[k]): $k of $(length(sacs_hp))")
	# τ_init_lp = 1.0/edgepitch(f0s[k], 0.060; kind=:lp)
	# τ_init_hp = 1.0/edgepitch(f0s[k], 0.060; kind=:hp)
	τ_init_lp = 1.0/f0s[k]
	τ_init_hp = 1.0/f0s[k]
	# est_lp = Adapt.est_τ0(τ, τ_init_lp, sacs_lp[k], sieve;
	# 					  μ_init=1e-8, verbose=true)
	# est_hp = Adapt.est_τ0(τ, τ_init_hp, sacs_hp[k], sieve;
	# 					  μ_init=1e-8, verbose=true)
	est_lp = Adapt.est_τ0(τ, τ_init_lp, a_lp(τ, f0s[k]), sieve;
						  μ_init=μ_est(f0s[k]), verbose=true)
	est_hp = Adapt.est_τ0(τ, τ_init_hp, a_hp(τ, f0s[k]), sieve;
						  μ_init=μ_est(f0s[k]), verbose=true)
	push!(ests_lp, est_lp)
	push!(ests_hp, est_hp)
end

p2 = plot(f0s, (1.0 ./ ests_hp ) ./ f0s, ylims=(0.9, 1.15))#, xscale=:ln)
plot!(f0s, (1.0 ./ ests_lp ) ./ f0s)

plot!(f0s, edgepitch.(f0s, 0.015; kind=:lp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.030; kind=:lp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.060; kind=:lp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)

plot!(f0s, edgepitch.(f0s, 0.015; kind=:hp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.030; kind=:hp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.060; kind=:hp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)

end
