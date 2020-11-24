module Salience

include("../util/PeakPick.jl")
include("../core/Sieves.jl")

using .PeakPick: findpeaknear, peakinterp

function salience(τ::Array{Float64,1},
				  τ0s::Array{Float64,1},
				  input,
				  sieve;
				  normalize=true)
	@assert length(τ) == length(input) "delays and input must be of same length for salience calc"
	sal = zeros(length(τ0s))
	for (k,τ0) ∈ enumerate(τ0s)
		if τ0 == 0.0
			@inbounds sal[k] = 0.0
			continue
		end
		s = sieve(τ, τ0)
		if normalize
			#norm = sqrt(sum( s .^ 2))
			norm = sum( abs.(s) ) 
			if norm == 0.0
				@inbounds sal[k] = 0.0
				continue
			end
			s ./= norm
		end
		@inbounds sal[k] = s' * input
	end
	return sal
end

"""
"""
function peak_sal(t, t0, x;
				  boundf=PeakPick.bounds_trough,
				  sievef= (t, t0) -> Sieves.sievegauss(t, t0, 1e-4))
    (i_lo, i_hi) = boundf(t, t0, x) 	# this should be reconsidered, 
										# since salience peak global property
	if i_hi - i_lo < 2
		i_lo -= 1
		i_hi += 1
	end
    sal = salience(t, t[i_lo:i_hi], x, sievef)
	i_max = argmax(sal)
	if i_max >= length(sal)
		i_max = length(sal) - 1
	end
	if i_max <= 1
		i_max = 2
	end
	i_maxo = i_max + i_lo - 1
	t_max = PeakPick.peakinterp(t[i_maxo-1:i_maxo+1], sal[i_max-1:i_max+1])
    return t_max
end


function est_τ0(τ, τ0s, input, sieve; τ0_init=:none, normalize=true)
	if τ0_init == :none
		τ0_init = 1.0 / ( (1/τ0s[end] + 1/τ0s[1]) /2 ) 
	end
	#s = salience( τ, τ0s, input, sieve; normalize=normalize)
	# idx = findpeaknear( τ0s, s, τ0_init; tol=0.25 )
	## idx = argmax(s)
	## if idx == 1
	## 	idx += 1
	## elseif idx == length(s)
	## 	idx -= 1
	## end
	## τ0_est = peakinterp( τ0s[idx-1:idx+1], s[idx-1:idx+1] )
	#PeakPick.peak_strict(τ0s, τ0_init, s)
	peak_sal(τ, τ0_init, input; sievef=sieve)
end

end
