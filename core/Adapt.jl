module Adapt

using Plots

"""
	d(τ, τ0, sieve; dτ=0.00001, norm=true)

Takes the numerical derivative of a sieve function with respect to τ0.  norm determines whether to rescale the resulting function so that its maximum is equal to 1.0
"""
function d(τ, τ0, sieve; dτ=0.00005, norm=true)
	sp = sieve(τ, τ0 + dτ)
	sn = sieve(τ, τ0 - dτ)
	ds = (sp .- sn) ./ (2.0*dτ)
	# if norm
	# 	maxval = maximum(abs.(ds))
	# 	if maxval > 0.0
	# 		ds ./= maxval 
	# 	end
	# end
	return ds
end


"""
	adapt(τ, τ0, input)

Given an input and an estimate fundamental delay value, it updates the estimate.
"""
function adapt(τ, τ0_init, input, sieve; μ=1e-5, maxchange = 0.1)
	τ0_new = τ0_init + μ * ( d(τ, τ0_init, sieve; dτ=μ)' * input )
end

## Old adaptive step size stuff
#function adapt(τ, τ0_init, input, sieve; μ=5e-13, maxchange = 1.1, tol=0.00001)
# 	τ0_0 = τ0_init + μ * ( d(τ, τ0_init, sieve; dτ=μ)' * input )
# 	τ0_1 = τ0_init + (μ/2.0) * ( d(τ, τ0_init, sieve; dτ=μ/2)' * input )
# 	τ0_1 = τ0_1 + (μ/2.0) * ( d(τ, τ0_1, sieve; dτ=μ/2)' * input )
# 	err = τ0_1 - τ0_0
# 	# println("τ0_0: ", τ0_0)
# 	# println("τ0_1: ", τ0_1)
# 	# println("err: ", err)
# 	μ_new = 0.9 * μ * min( max( tol / abs(err), 0.3), 2.0)
# 	# println("New μ: ", μ_new)
# 	if τ0_1 <= 0.0 error("τ0 too small!") end
# 	return (τ0_1, μ_new)
# end
"""
	est_τ0(τ, τ0_init, input, sieve; μ=1e-5, ϵ=1e-9, maxits=100, anim=false, verbose=false)

Estimate the fundamental delay. sieve must be a function taking τ and τ0.
"""
function est_τ0(τ, τ0_init, input, sieve; μ_init=5e-13, ϵ=1e-9, maxits=100, 
				anim=false, verbose=false, xlims=(0.0, 0.03))
	τ0_est = τ0_init
	τ0_new = 0.0
	converged = false
	count = 0
	μ = μ_init
	if anim
		@gif for k ∈ 1:maxits
			count += 1
			#(τ0_new, μ) = adapt(τ, τ0_est, input, sieve; μ=μ)
			τ0_new = adapt(τ, τ0_est, input, sieve; μ=μ)
			if abs(τ0_new - τ0_est) < ϵ
				converged = true
				break
			end
			τ0_est = τ0_new
			plot(τ, input, xlims=xlims)
			plot!(τ, sieve(τ, τ0_est))
		end
		if verbose println("Converged after $count steps: $converged") end
	else
		for k ∈ 1:maxits
			count += 1
			# (τ0_new, μ) = adapt(τ, τ0_est, input, sieve; μ=μ)
			τ0_new = adapt(τ, τ0_est, input, sieve; μ=μ)
			if abs(τ0_new - τ0_est) < ϵ
				converged = true
				break
			end
			τ0_est = τ0_new
		end
		if verbose println("Converged after $count steps: $converged") end
	end
	return τ0_new
end

end
