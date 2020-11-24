module EdgePitchUtils

include("./PeakPick.jl")
using .PeakPick: peakinterp, findpeaknear

################################################################################
# Edge pitch theory (from Hartmann et al.)
################################################################################
"""
	edgepitch(fe, τ_max; kind=:lp)

Get theoretical value of edge pitch. τ_max is the maximum delay.  Available kinds are :hp and :lp.
"""
function edgepitch(fe, τ_max; kind=:lp)
	n = Int(round(fe*τ_max))
	γ = 0.57722 # Euler-Mascheroni constant
	if kind == :lp
		p = fe / ( 1 + (1 / (4 * n)) * (log(n) + γ + 1 / (2 * n)) )
	elseif kind == :hp
		p = fe / ( 1 - (1 / (4 * n)) * (log(n) + γ + 1 / (2 * n)) )
	end
end

"""
	a_lp(τ, fe)

Sinc approximation of SAC for low-pass edge pitch stimulus
"""
a_lp(τ, fe) = sinc.(2*fe .* τ)


"""
	a_hp(τ, fe)

Sinc approximation of SAC for high-pass edge pitch stimulus
"""
a_hp(τ, fe) = 1.0 .- sinc.(2*fe .* τ)


"""
	delayfrompeaks(peaks; maxpeaks=4)

Calculate the fundamental delay time from a series of peaks (all equally weighted).
"""
function delayfrompeaks(peaks; maxpeaks=4)
	N = min(length(peaks), maxpeaks)
	τ_est = (1.0/N) * sum( [ peaks[n]/n for n ∈ 1:min(length(peaks), maxpeaks) ] )
end


"""
	delayfrompeaks(peaks; maxpeaks=4)

Calculate the fundamental delay time from a series of peaks (all equally weighted).
"""
function delayfrompeaks(peaks; maxpeaks=4, envf=x->1.0)
	N = min(length(peaks), maxpeaks)
	coefs = [envf(p) for p ∈ peaks[1:N]]
	norm = sum(coefs)
	τ_est = (1.0/norm) * sum( [coefs[n]*peaks[n]/n for n ∈ 1:N] )
end



"""

Calculate the fundamental delay time from a series of peaks (weighted by peak height).
"""
function delayfrompeaks(peaks, peakvals; maxpeaks=4)
	N = min(length(peaks), maxpeaks)
	τ_est = (1.0/N) * sum( [ peakvals[n]*peaks[n]/n for n ∈ 1:min(length(peaks), maxpeaks) ] )
end


end
