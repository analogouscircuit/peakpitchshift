module Sieves

using Match
include("Bumps.jl")
using .Bumps: gauss, powgauss, hat, laplace

export sievegauss, sievehat

"""
	sievegauss(τ, τ0, σ; skew=0.0, env=:none)

Gaussian sieve.  σ should be an array, containing the width of each bump.  The number of values determines the number of bumps. The env (envelope) option takes a function that itself takes two arguments: (τ, τ0).
"""
function sievepowgauss(τ, τ0, σ, γ=0.5; skew=0.0, env=:none)
	sieve = sum([powgauss(τ, p*τ0, σ[p], γ; skew=skew) for p in 1:length(σ)])
	if env != :none
		sieve .*= env(τ, τ0)
	end
	return sieve
end


"""
	sievegauss(τ, τ0, σ; skew=0.0, env=:none)

Gaussian sieve.  σ should be an array, containing the width of each bump.  The number of values determines the number of bumps. The env (envelope) option takes a function that itself takes two arguments: (τ, τ0).
"""
function sievegauss(τ, τ0, σ; skew=0.0, env=:none)
	sieve = sum([gauss(τ, p*τ0, σ[p]; skew=skew) for p in 1:length(σ)])
	if env != :none
		sieve .*= env(τ, τ0)
	end
	return sieve
end


"""
	sievefunc(τ, τ0, σ; bump=:gauss skew=0.0, env=:none)
"""
function sievefunc(τ, τ0, σ; bump=:gauss, skew=0.0, env=:none)
	bumpf(t, t0, s) = @match bump begin
		:gauss 		=> Bumps.gauss(t, t0, s; skew=skew)
		:hat 		=> Bumps.hat(t, t0, s; skew=skew)
		:powgauss 	=> Bumps.powgauss(t, t0, s, 0.5; skew=skew)
		:laplace 	=> Bumps.laplace(t, t0, s)
		:cosbump 	=> Bumps.cosbump(t, t0, s)
		:cosbump2 	=> Bumps.cosbump2(t, t0, s)
		:cosbump4 	=> Bumps.cosbump4(t, t0, s)
		:cosbump6 	=> Bumps.cosbump6(t, t0, s)
		:cosbump0 	=> Bumps.cosbump0(t, t0, s)
		:rect 		=> Bumps.rect(t, t0, s)
		:gauss2 	=> Bumps.powgauss(t, t0, s, 2.0)
		:gauss4 	=> Bumps.powgauss(t, t0, s, 4.0)
	end
	sieve = sum([bumpf(τ, p*τ0, σ[p]) for p in 1:length(σ)])
	if env != :none
		sieve .*= env(τ, τ0)
	end
	return sieve
end


"""
	sievelaplace(τ, τ0, σ; skew=0.0, env=:none)

Gaussian sieve.  σ should be an array, containing the width of each bump.  The number of values determines the number of bumps. The env (envelope) option takes a function that itself takes two arguments: (τ, τ0).
"""
function sievelaplace(τ, τ0, σ; skew=0.0, env=:none)
	sieve = sum([laplace(τ, p*τ0, σ[p]) for p in 1:length(σ)])
	if env != :none
		sieve .*= env(τ, τ0)
	end
	return sieve
end

"""
	sievehat(τ, τ0, σ; skew=0.0, dτ=0.00001, env=:none)

Mexican hat sieve.  σ can be an array.  The number of bumps generated will correspond to the length of σ. The envelope option takes a function which should itself take two arguments: (τ, τ0). 
"""
function sievehat(τ, τ0, σ; skew=0.0, env=:none)
	sieve = sum([hat(τ, p*τ0, σ[p]; skew=skew) for p in 1:length(σ)])
	if env != :none
		sieve .*= env(τ, τ0)
	end
	return sieve
end

end
