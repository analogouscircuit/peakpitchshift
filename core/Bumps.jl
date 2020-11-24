module Bumps

using SpecialFunctions: erf

export gauss, hat


"""
	skewfunc(x, μ, σ, skew)

Skew function.  Multiple this by a Gaussian scaled up by a factor of two (but with the same parameters) to produce a skew Gaussian. Skew 
"""
function skewfunc(x, μ, σ, skew)
	0.5 .* (1.0 .+ erf.( (skew/√2) .* ((x .- μ) ./ σ ) ))
end


"""
"""
function laplace(τ, τ0, σ)
	exp.( -abs.( τ .- τ0 ) ./ σ)
end

"""
	gauss(τ, τ0, σ; skew=0.0)

Produces an unnormalized Gaussian bump (height of mean is 1). Option to skew.
"""
function gauss(τ, τ0, σ; skew=0.0)
	arg  = -0.5 .* ((τ .- τ0) .^ 2) ./  (σ)^2
	bump = exp.( arg )
	if skew != 0.0
		bump .*= 2 .* skewfunc(τ, τ0, σ, skew) 
		bump ./= maximum(bump)
	end
	return bump
end

"""
	powgauss(τ, τ0, σ, γ; skew=0.0)

Like gauss, except that the resulting bump is raised to a power (and thus sharpened or broadened).
"""
function powgauss(τ, τ0, σ, γ=0.5; skew=0.0)
	arg  = -0.5 .* ((τ .- τ0) .^ 2) ./  (σ)^2
	bump = exp.( arg )
	if skew != 0.0
		bump .*= 2 .* skewfunc(τ, τ0, σ, skew) 
		# bump ./= maximum(bump)
	end
	return bump .^ γ
end


"""
	hat(τ, τ0, σ; skew=0.0, dτ=0.000001)

Produces a "mexican hat" function (or Ricker wavelet).  Does so by calculating the second derivative of a Gaussian numerically.  Numerical approach is chosen to allow for skew (otherwise could of course use analytical expression). 
"""
# Analytic expression (not strictly negative second derivative of skewed
# Gaussian)
function hat(τ, τ0, σ; skew=0.0)
	h = ( 2 / (√(3σ) * π^(1/4)) ) .* 
		( 1.0 .- ((τ .- τ0) ./ σ) .^ 2) .* 
		exp.( - ((τ .- τ0) .^ 2) ./ (2 * σ ^2) )
	if skew != 0.0
		h .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	h ./= maximum(h)
	return h
end

function cosbump(t, t0, σ; skew=0.0)
	args = (2π/σ) .* (t .- t0)
	out = map( x -> x < -π/2 || x > π/2 ? π/2 : x, args ) .|> cos
	if skew != 0.0
		out .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	out
end

function cosbump2(t, t0, σ; skew=0.0)
	args = (2π/σ) .* (t .- t0)
	out = map( x -> x < -π/2 || x > π/2 ? π/2 : x, args ) .|> cos .|> x -> x^2
	if skew != 0.0
		out .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	out
end

function cosbump4(t, t0, σ; skew=0.0)
	args = (2π/σ) .* (t .- t0)
	out = map( x -> x < -π/2 || x > π/2 ? π/2 : x, args ) .|> cos .|> x -> x^4
	if skew != 0.0
		out .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	out
end

function cosbump6(t, t0, σ; skew=0.0)
	args = (2π/σ) .* (t .- t0)
	out = map( x -> x < -π/2 || x > π/2 ? π/2 : x, args ) .|> cos .|> x -> x^6
	if skew != 0.0
		out .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	out
end

function cosbump0(t, t0, σ; skew=0.0)
	args = (2π/σ) .* (t .- t0)
	out = map( x -> x < -π/2 || x > π/2 ? π/2 : x, args ) .|> cos .|> x -> x^0.5
	if skew != 0.0
		out .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	out
end

function rect(t, t0, σ; skew=0.0)
	rectfunc(τ, w) = abs(τ) < w ? 1.0 : 0.0		
	out = map(x -> rectfunc(x-t0, σ), t)
	if skew != 0.0
		out .*= 2 .* skewfunc(τ, τ0, σ, skew)
	end
	out
end


end
