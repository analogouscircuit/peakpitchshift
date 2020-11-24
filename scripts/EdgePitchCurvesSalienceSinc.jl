module EdgePitchCurves

using JLD
using Plots; pyplot()

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
fs = 96000
τ = collect(range(0, stop=0.1, step=1/fs))

σ0 = 5e-4 #width at τ0 = 1/155
sig(τ0) = Sigmas.uniform_scaled(τ0; σ0=σ0, num_bumps=60)
#sig(τ0) = Sigmas.nonuniform_scaled(τ0; σ0=σ0, num_bumps=40)
#sig(τ0) = Sigmas.uniform(; σ0=σ0, num_bumps=30)
# env = Envelopes.hcc
#env = Envelopes.gamma
#env(τ, τ0) = Envelopes.gamma(τ, τ0*5.0; peakstart=true)
env = Envelopes.oneovern
#env = Envelopes.tempift
#sig(τ0) = Sigmas.nonuniform_envinv(σ0, env; σ0=0.0, num_bumps=5)
sieve(τ, τ0) = Sieves.sievegauss(τ, τ0, sig(τ0); env=env)

ests_lp = []
ests_hp = []
for k in 1:length(f0s)
	println("$(f0s[k]): $k of $(length(f0s))")
	#τ_init_lp = 1.0/edgepitch(f0s[k], 0.060; kind=:lp)
	#τ_init_hp = 1.0/edgepitch(f0s[k], 0.060; kind=:hp)
	τ_init_lp = 1.0/f0s[k]
	τ_init_hp = 1.0/f0s[k]
	#τ0s = collect(range(0.75/f0s[k], stop=1.25/f0s[k], step=1/fs))
	#τ0s = filter( x -> x >= 0.75/f0s[k] && x <= 1.25/f0s[k], τ)  
	# s_lp = Salience.salience(τ, τ0s, a_lp(τ, f0s[k]), sieve)
	# s_hp = Salience.salience(τ, τ0s, a_hp(τ, f0s[k]), sieve)
	# est_lp = Salience.est_τ0(τ, τ0s, a_lp(τ, f0s[k]), sieve; τ0_init = τ_init_lp, normalize=true)
	# est_hp = Salience.est_τ0(τ, τ0s, a_hp(τ, f0s[k]), sieve; τ0_init = τ_init_hp, normalize=true)
	est_lp = Salience.peak_sal(τ, τ_init_lp, a_lp(τ, f0s[k]); sievef=sieve)
	est_hp = Salience.peak_sal(τ, τ_init_hp, a_hp(τ, f0s[k]); sievef=sieve)
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
