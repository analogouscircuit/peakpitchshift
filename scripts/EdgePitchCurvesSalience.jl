module EdgePitchCurves

using JLD
using Plots

include("../core/Envelopes.jl")
include("../core/Sigmas.jl")
include("../core/Sieves.jl")
include("../core/Salience.jl")
include("../util/EdgePitchUtils.jl")

using .EdgePitchUtils: edgepitch, a_lp, a_hp

################################################################################
#
################################################################################


################################################################################
# Main Script
################################################################################
data = load("/home/dahlbom/research/AdaptiveSieves/stimuli/sac_edge_data_avg.jld")
f0s = data["f0s"]
fs = data["fs"]
sacs_lp = data["sacs_lp"]
sacs_hp = data["sacs_hp"]
τ = data["τ"] 

σ0 = 5e-4 #width at τ0 = 1/155
sig(τ0) = Sigmas.uniform_scaled(τ0; σ0=σ0, num_bumps=40)
#sig(τ0) = Sigmas.nonuniform_scaled(τ0; σ0=σ0, num_bumps=40)
#sig(τ0) = Sigmas.uniform(; σ0=σ0, num_bumps=30)
env = Envelopes.hcc
# env = Envelopes.gamma
#env = Envelopes.tempift
#sig(τ0) = Sigmas.nonuniform_envinv(σ0, env; σ0=0.0, num_bumps=5)
sieve(τ, τ0) = Sieves.sievegauss(τ, τ0, sig(τ0); env=env)

ests_lp = []
ests_hp = []
for k in 1:length(f0s)
	println("$(f0s[k]): $k of $(length(f0s))")
	# τ_init_lp = 1.0/edgepitch(f0s[k], 0.060; kind=:lp)
	# τ_init_hp = 1.0/edgepitch(f0s[k], 0.060; kind=:hp)
	τ_init_lp = 1.0/f0s[k]
	τ_init_hp = 1.0/f0s[k]
	τ0s = collect(range(0.75/f0s[k], stop=1.25/f0s[k], step=1/fs))
	est_lp = Salience.est_τ0(τ, τ0s, sacs_lp[k], sieve; τ0_init = τ_init_lp)
	est_hp = Salience.est_τ0(τ, τ0s, sacs_hp[k], sieve; τ0_init = τ_init_hp)
	push!(ests_lp, est_lp)
	push!(ests_hp, est_hp)
end

l = @layout [ a{0.8h}; b]

p1 = plot(f0s, (1.0 ./ ests_hp ) ./ f0s, ylims=(0.8, 1.2))#, xscale=:ln)
plot!(f0s, (1.0 ./ ests_lp ) ./ f0s)

plot!(f0s, edgepitch.(f0s, 0.015; kind=:lp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.030; kind=:lp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.060; kind=:lp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)

plot!(f0s, edgepitch.(f0s, 0.015; kind=:hp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.030; kind=:hp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)
plot!(f0s, edgepitch.(f0s, 0.060; kind=:hp) ./ f0s, linestyle=:dash, linecolor=:black, linewidth=.5)

p2 = plot(τ, sieve(τ, 1/f0s[1]))
plot!(τ, sieve(τ, 1/f0s[25]))

p = plot(p1, p2, layout=l)

end
