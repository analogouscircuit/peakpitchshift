module Envelopes

export gamma, tempift

"""
	hcc(t, τ0)

	Hartmann-Cariani-Colburn (2019) envelope.  All values are 1.0 until cutoff
"""
function hcc(t, τ0)
	τ200 = 0.0556
	τ2500 = 0.0139
	τ_max(f) = τ2500 + (τ200 - τ2500) * log(2500.0/f) / log(2500.0/200)
	fe = max(1e-32, 1/τ0)
	t_out = map( x -> x < τ_max(fe) ? 1.0 : 0.0, t)
end

function hcc_maxdelay(τ0; τ200=0.0556, τ2500=0.0139)
	# τ200 = 0.0556
	# τ2500 = 0.0139
	τ_max(f) = τ2500 + (τ200 - τ2500) * log(2500.0/f) / log(2500.0/200)
	fe = max(1e-32, 1/τ0)
	τ_max(fe)
end

"""
	gamma(t, τ0; n=4)

Envelope of a gamma chirp. By default, the envelope starts at the peak of the chirp. n is the order of the filter.
"""
function gamma(t, τ0; n=4, peakstart = false)
	f0 = 1.0/τ0
	B = 24.7 + 0.108*f0
	b = 1.018 * B
	if peakstart
		t0 = (n-1)/(2π * b) # take only from peak
	else
		t0 = 0.0
	end
	env = (t .+ t0) .^ (n-1) .* exp.(-2π*b .* (t .+ t0))
	# env ./= maximum(env)
end

function oneovern(t, τ0)
	env = τ0 ./ t
	map( x -> x == NaN || x == Inf ? 0.0 : x, env)
end

"""
	expenv(t, τ_decay)

Exponential envelope
"""
function expenv(t, τ_decay)
	exp.(-t ./ τ_decay)
end

function oneoverrootn(t, τ0)
	env = sqrt.(τ0 ./ t)
	map( x -> x == NaN || x == Inf ? 0.0 : x, env)
end

"""
	tempift(τ, τ0; num_h=5, α=0.031)

Envelope function derived from taking the inverse Fourier transform of a harmonic template in spectral domain.
"""
function tempift(τ, τ0; num_h=5, α=0.031)
	env = sum([(2π*n*α/τ0) .* exp.( -2.0 * (π*n*α/τ0)^2 .* (τ .^ 2) ) for n in 1:num_h]) 
	#env ./= maximum(env)
end


"""
	w_berox(l, cf; nc=10.8, nd=2.0, a=200.0)

Envelope from Bernstein and Oxenham. Only good down to 1500 Hz.
"""
function w_berox(l, cf; nc=10.8, nd=2.0, a=200.0)
	# Bernstein and Oxenham
	# only good for freqs > 1500.0
	cf0 = 1500.0 # Hz
	l0 = 0.033 #seconds
	m = (cf^2 - a + (a/l0) * ( (nc+nd) / cf) ) / (nd / cf)
	b1 = 0.5/cf
	b2 = nc/cf
	b3 = (nc + nd)/cf
	println("Boundaries: $b1 $b2 $b3")
	if l < b1
		return 0.0
	elseif (b1 <= l) && (l < b2)
		return (cf ^ 2.0) / cf0
	elseif (b2 <= l) && (l < b3)
		return (cf^2 / cf0) - m * (l - nc/cf)
	elseif b3 <= l
		return a - (a/l0) * l
	else
		error("something logically impossible has occurred...")
	end
end


"""
	w_pm(t; w =-69.0, tp=0.0055, ts=1.9)

Envelope from Plack and Moore.  Double exponential.  Sample parameter values available in Envelopes.jl.
"""
function w_pm(t; w =-69.0, tp=0.0055, ts=1.9)
	#convert from dB to scalar
	w = 10 ^ (w/10)
	w_f(τ) = (1.0 - w) * ( 1 + 2τ / tp) * exp( - 2τ / tp ) +
			 w * (1.0 + 2τ / ts) * exp( - 2τ/ ts)
	w_f.(t)
end

## Sample parameter values
freqs = [300.0, 900.0, 2700.0, 8100.0]
dBs = [44.0 54.0 64.0; 32.0 42.0 52.0; 20.0 30.0 40.0; 29.0 39.0 49.0]
# idx1: freq, idx2: level (as referenced above)
t_pb = [8.6 7.6 7.1; 6.8 6.0 5.1; 6.3 5.4 4.7; 5.9 4.9 3.8] .* 0.001
t_pa = [4.5 4.9 5.5; 4.0 3.5 2.5; 2.9 2.6 2.4; 2.9 2.9 2.1] .* 0.001
w    = [-50.0 -56.0 -69.0; -39.0 -51.0 -51.0; -34.0 -39.0 -44.0; -33.0 -40.0 -46.0]
t_sb = [4700.0 5200.0 3100.0; 48.0 57.0 27.0; 36.0 30.0 27.0; 36.0 33.0 26.0] .* 0.001 
t_sa = [2700.0 3600.0 1900.0; 28.0 32.0 13.0; 17.0 14.0 14.0; 17.0 18.0 14.0] .* 0.001

end
