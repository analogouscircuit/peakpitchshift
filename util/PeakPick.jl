module PeakPick


include("../core/Adapt.jl")
include("../core/Sieves.jl")
#include("../core/Salience.jl")


export interpfunc, quadinterp, peakinterp, grad1D, crossings, findpeaknear, sparsify

function interpfunc(x, y)
	if length(x) != length(y) error("Different number of x and y points!") end
	ord = length(x) - 1		
	l_funcs = 
		[l(z) = prod( [(z - x[m])/(x[j] - x[m]) for m ∈ 1:length(x) if m != j] ) for j ∈ 1:length(x) ]
	f(z) = sum( [y[j]*l_funcs[j](z) for j ∈ 1:length(x)] )
end

function quadinterp(x, y)
      if length(x) != 3 && length(y) != 3
            error("Need three points")
      end
      l1(x_in) = (x_in - x[2]) * (x_in - x[3]) / ( (x[1] - x[2]) * (x[1] - x[3]) )
      l2(x_in) = (x_in - x[1]) * (x_in - x[3]) / ( (x[2] - x[1]) * (x[2] - x[3]) )
      l3(x_in) = (x_in - x[1]) * (x_in - x[2]) / ( (x[3] - x[1]) * (x[3] - x[2]) )
      f(x_in) = y[1]*l1(x_in) + y[2]*l2(x_in) + y[3]*l3(x_in)
      return f
end

function peakinterp(x, y; testpeak=false)
      if length(x) != 3 && length(y) != 3
            error("Need three points")
      end
	  # if three points don't surround a peak, just send back maximum
	  if testpeak
		  if (y[1] < y[2] < y[3]) || (y[1] > y[2] > y[3])
			  idx = argmax(y)
			  return x[idx]
		  end
	  end
      k1 = y[1] / ( (x[1] - x[2]) * (x[1] - x[3]) )
      k2 = y[2] / ( (x[2] - x[1]) * (x[2] - x[3]) )
      k3 = y[3] / ( (x[3] - x[1]) * (x[3] - x[2]) )
      return (k1 * (x[2] + x[3]) + k2 * (x[1] + x[3]) + k3 * (x[1] + x[2])) /
                  (2 * (k1 + k2 + k3) )
end


function grad1D(x, h)
	diff3p(x, h) = (x[3] - x[1]) / (2h)
	diff2p(x, h) = (x[2] - x[1]) / h
      len = length(x)
      out = Array{typeof(x[1]), 1}(undef, len)
      out[1] = diff2p(x[1:2], h)
      for k ∈ 2:(len-1)
            out[k] = diff3p(x[k-1:k+1], h)
      end
      out[end] = diff2p(x[end-1:end], h)
      return out
end


function crossings(x; kind=:max)
	idcs = []
	if kind == :max
		for k ∈ 1:length(x)-1
			if x[k] >= 0 && x[k+1] < 0
				  push!(idcs, k)
			end
		end
	elseif kind == :both
		for k ∈ 1:length(x)-1
			if (x[k] >= 0 && x[k+1] < 0) || (x[k] < 0 && x[k+1] >=0)
				  push!(idcs, k)
			end
		end
	elseif kind == :min
		for k ∈ 1:length(x)-1
			if x[k] < 0 && x[k+1] >=0
				  push!(idcs, k)
			end
		end
	else
		error("kind must be :max, :min or :both")
	end
	return idcs
end


function findpeaknear(x, y, val; tol=0.025)
	dx = x[2] - x[1] 	# assumes equal spacing on x-axis
	grad = PeakPick.grad1D(y, dx)
	idcs = PeakPick.crossings(grad)
	idx = 0
	for k ∈ idcs
		if abs(1.0 - x[k]/val) < tol 
			idx = k
			break
		end
	end
	if idx == 0
		println("No zero-crossing near expected value")
		return 0
	end
	return idx
end

function sparsify(τ, sac)
	crossidcs = PeakPick.crossings(PeakPick.grad1D(sac, τ[2]-τ[1]); kind=:max)
	τ_vals = zeros(length(crossidcs))
	s_vals = zeros(length(crossidcs))
	for (k, idx) in enumerate(crossidcs)
		τ_vals[k] = PeakPick.peakinterp( τ[idx-1:idx+1], sac[idx-1:idx+1] )
		f = PeakPick.quadinterp(τ[idx-1:idx+1], sac[idx-1:idx+1])
		s_vals[k] = f(τ_vals[k])
	end
	return (τ_vals, s_vals) 
end


# """
# 	hcpeaks(sac, τ0)
# 
# Find the peaks near subharmonics of τ0.  
# """
# function hcpeaks(τ, τ0, sac; tol=0.1)
# 	dt = τ[2] - τ[1]
# 	fs = 1/dt
# 	peaks = []
# 	k = 1
# 	while k * τ0 <= (τ[end] - τ0/2)
# 		# idx = findpeaknear(τ, sac, k*τ0; tol=tol) 
# 		τ_min = k*τ0*(1.0 - tol)
# 		τ_max = k*τ0*(1.0 + tol)
# 		imin = max( 1, Int(floor(τ_min*fs)) )
# 		imax = min( length(sac), Int(ceil(τ_max*fs)) )
# 		idx = argmax(sac[imin:imax])
# 		idx += imin - 1
# 		if idx <= 1 idx = 2 end
# 		if idx >= length(sac) idx = length(sac) - 1 end
# 		push!(peaks, peakinterp(τ[idx-1:idx+1], sac[idx-1:idx+1]))
# 		k += 1
# 	end
# 	return peaks
# end


"""
	bounds_infl(t, x, t0)

Give the indices of the nearest inflection points on either side of a starting point (t0).  The second derivative of x(t0) must be negative, i.e. it must be near a maximum.
"""
function bounds_infl(t, t0, x)
    @assert length(t) == length(x)
    n = length(t)
    dt = t[2] - t[1]
    fs = 1/dt
	idx = Int(round((t0-t[1])*fs))
    d1 = PeakPick.grad1D(x, dt)
    d2 = PeakPick.grad1D(d1, dt)
    if d2[idx] >= 0 error("Not near a maximum!") end
    k = 0
    i_hi = 0
    i_lo = 0
    while idx + k < n
         if d2[idx+k] <= 0.0 && d2[idx+k+1] > 0.0
            i_hi = idx + k
            break
        end
        k += 1
    end
    if i_hi == 0; i_hi = n end
    k = 0
    while idx - k > 1
         if d2[idx-k] <= 0 && d2[idx-k-1] > 0
            i_lo = idx - k
            break
        end
        k += 1
    end
    if i_lo == 0; i_lo = 1 end
    return (i_lo, i_hi) 
end


"""
	bounds_trough(t, x, t0)

Give the indices of the nearest minima on either side of a starting point (t0).  The second derivative of x must be negative near t0, i.e. x(t0) must be near a maximum.
"""
function bounds_trough(t, t0, x)
    @assert length(t) == length(x) "Input and times must be the same length"
	@assert (t[1] <= t0) && (t[end] >= t0) "τ0 not within given time series" 
    n = length(t)
    dt = t[2] - t[1]
    fs = 1/dt
	idx = Int(round(t0*fs - t[1]*fs))
    d1 = PeakPick.grad1D(x, dt)
    # d2 = PeakPick.grad1D(d1[idx-3:idx+3], dt)
    # if d2[4] >= 0 error("Not near a maximum!") end
    k = 0
    i_hi = 0
    i_lo = 0
    while idx + k < n
         if d1[idx+k] <= 0.0 && d1[idx+k+1] > 0.0
            i_hi = idx + k
            break
        end
        k += 1
    end
    if i_hi == 0; i_hi = n end
    k = 0
    while idx - k > 1
         if d1[idx-k] >= 0 && d1[idx-k-1] < 0
            i_lo = idx - k
            break
        end
        k += 1
    end
    if i_lo == 0; i_lo = 1 end
    return (i_lo, i_hi) 
end

"""
"""
function getdelays(τ0, τ_max)
    n_max = Int(floor( (τ_max - τ0/2) / τ0))
    [τ0*k for k in 1:n_max]
end

"""
"""
function peak_strict(t, t0, x; boundf=PeakPick.bounds_trough)
    (i_lo, i_hi) = boundf(t, t0, x)
    i_max = argmax(x[i_lo:i_hi]) + i_lo - 1
	i_l = max(1, i_max - 1)
	i_h = min(length(t), i_max + 1)
	if i_h - i_l > 1
		t_max = PeakPick.peakinterp(t[i_l:i_h], x[i_l:i_h]) 
	elseif i_h - i_l == 1
		t_max = (x[i_l]*t[i_l] + x[i_h]*t[i_h]) / (x[i_l] + x[i_h])
	else
		t_max = t[i_max]
	end
    return t_max
end

"""
"""
function peak_int(t, t0, x; boundf=PeakPick.bounds_trough)
    m = minimum(x)
    if m < 0.0; x .+= abs(m) end
    (i_lo, i_hi) = boundf(t, t0, x)
    t_max = (t[i_lo:i_hi]' * x[i_lo:i_hi]) / sum(x[i_lo:i_hi])
    return t_max
end


"""
"""
function peak_gauss(t, t0, x; 
					boundf=PeakPick.bounds_trough, σ=4e-5, μ=1e-8, 
					verbose=false, anim=false, xlims=(0, 0.01))
    μ /= σ
    (i_lo, i_hi) = boundf(t, t0, x)
    i_max = argmax(x[i_lo:i_hi]) + i_lo - 1
    t0_init = t[i_max]
    sieve(τ, τ0) = Sieves.sievegauss(τ, τ0, σ; skew=0.0, env=:none)
    t0_est = Adapt.est_τ0(t, t0_init, x, sieve; μ_init=μ, verbose=verbose, anim=anim, xlims=xlims)
    return t0_est
end


"""
"""
# function peak_sal(t, t0, x; boundf=PeakPick.bounds_trough, σ=1e-4)
#     (i_lo, i_hi) = boundf(t, t0, x)
#     sieve(τ, τ0) = Sieves.sievegauss(τ, τ0, σ; skew=0.0, env=:none)
#     sal = Salience.salience(t, t[i_lo:i_hi], x, sieve)
#     i_max = argmax(sal) + i_lo - 1
#     t_max = PeakPick.peakinterp(t[i_max-1:i_max+1], x[i_max-1:i_max+1])
#     return t_max
# end



"""
"""
function hcpeaks(t, t0, x; peakf = peak_strict)
    delays = getdelays(t0, t[end])
    [peakf(t, τ, x) for τ ∈ delays]
end
    

"""
"""
function peakdeviations(τ0, peaks)
    delays = [τ0*k for k ∈ 1:length(peaks)]
    return ( (1.0 ./ peaks) ./ (1.0 ./ delays) ) .- 1.0
end


"""
"""
function peakdeviations_hc(τ0, peaks)
    delays = [τ0*k for k ∈ 1:length(peaks)]
	peaks_rs = [peak/k for (k, peak) ∈ enumerate(peaks)]
    return ( (1.0 ./ peaks_rs) ./ (1.0 ./ τ0) ) .- 1.0
end


end
