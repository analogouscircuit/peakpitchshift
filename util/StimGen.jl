module StimGen

using AuditoryFilters
using DSP
using FFTW
using Distributed
using Random

include("ACUtils.jl")
using .ACUtils

"""
	energy(sig)

Calculate signal energy (use for normalization).
"""
function energy(sig)
	sum( sig' * sig) / length(sig)
end


"""
	hc(t, f0; amps=[1.0 for k ∈ 1:12], 

Generate a harmonic complex.
"""
function hc(dur, f0, fs; amps=[1.0 for k ∈ 1:12], 
						 mistunings=[1.0 for k ∈ 1:12],
						 phase=[0.0 for k ∈ 1:12])
	t = collect(range(0, stop=dur, step=1/fs))
	sig = sum([amps[p] .* cos.( 2π*f0*mistunings[p]*p .* t .+ phase[p]) for p ∈ 1:length(amps)])
	sig ./= maximum(abs.(sig))
end


"""
	hc(t, f0; amps=[1.0 for k ∈ 1:12], 

Generate a harmonic complex.
"""
function hcfm(dur, f0, fs, ffm, fp; amps=[1.0 for k ∈ 1:12], 
						            mistunings=[1.0 for k ∈ 1:12],
						            phase=[0.0 for k ∈ 1:12])
	t = collect(range(0, stop=dur, step=1/fs))
	sig = sum([amps[p] .* cos.( 2π*f0*mistunings[p]*p .* t .+ phase[p]
							   .+ fp * f0 .* cos.(2*pi*ffm .* t)) 
			   for p ∈ 1:length(amps)])
	sig ./= maximum(abs.(sig))
end

"""
	bumptrain(f0, dur, fs)

Generate a periodic train of delta functions.  Places delta function in closest bin to ideal time.
"""
function bumptrain(f0, dur, fs)
    dt = 1/fs
    num_pts = Int(floor(dur/dt))
    t = collect(1:num_pts) .* dt
    out = zeros(length(t))
    k = 1
    while k/f0 <= dur
        idx = Int(round(fs*k/f0))
        out[idx] = 1.0
        k += 1
    end
    out
end

"""
	bumptrain_interp(f0, dur, fs)

Generate a periodic train of delta functions.  Distributes delta function between two bins nearest to ideal time.
"""
function bumptrain_interp(f0, dur, fs)
    dt = 1/fs
    num_pts = Int(floor(dur/dt)) 
    t = collect(1:num_pts) ./ fs
    out = zeros(length(t))
    k = 1
    while k/f0 <= dur 
        time = k/f0
        i_lo = Int(floor(fs*k/f0))
        if i_lo < length(t)
            w1 = (t[i_lo+1] - time) / dt 
            out[i_lo] = w1
            out[i_lo+1] = 1-w1
        else
            w1 = (t[i_lo] + dt - time) / dt 
            out[i_lo] = w1
        end
        k += 1
    end
    return out
end

"""
	bumphc(f0, dur, fs; num_h = 0, amps=:none, mt=:none

Generate a train of superimposed delta functions, the frequencies of which form a harmonioc complex (or harmonic complex with mistunings, if specified). 
"""
function bumphc(f0, dur, fs; num_h = 0, amps=:none, mt=:none)
    if num_h == 0
        num_h = length(amps)
    end
    if amps == :none
        amps = [1.0 for k ∈ 1:num_h]
    end
    if mt == :none
        mt = [1.0 for k ∈ 1:num_h]
    end
    dt = 1/fs
    n = Int(floor(dur/dt))
    out = zeros(n)
    for k ∈ 1:num_h
        out .+= amps[k] .* bumptrain(k*mt[k]*f0, dur, fs)
    end
    return out
end


"""
	hwr(x)

Half-wave rectify.
"""
hwr(x) = map(y -> y < 0 ? 0 : y, x)


"""
	sac(sig, fb, fs; actype=:fft, dur=0.0, norm=true)

Generate a summary autocorrelation (SAC).  Needs to be provided with a filterbank,
which can be generating using make_erb_filterbank from AuditoryFilters.
"""
function sac(sig, fb, fs; actype=:fft, dur=0.0, norm=true, compress = x -> x,
			 channels=false,
			 acchannels=false)
	nchan = length(fb.filters)
	sig_filt = filt(fb, sig)
	if dur != 0.0
		# take after transient if duration is less than signal input
		len_n = Int(round(dur*fs))
		len_sig = length(sig_filt[:,1])
		sig_filt = sig_filt[end-len_n+1:end,:]
	end
	sig_hwr = hwr(sig_filt)
	sig_hwr_ac = copy(sig_hwr)
	ac = actype == :fft ? ACUtils.ac_fft : ACUtils.ac_true
	for k ∈ 1:nchan
		sig_hwr_ac[:,k] = ac(sig_hwr_ac[:,k]) 
		sig_hwr_ac[:,k] = compress.(sig_hwr_ac[:,k])
	end
	out = sum(sig_hwr_ac, dims=2)[:,1]
	if norm
		#out ./= sqrt(energy(out))
		out ./= out[1] 
	end
	if channels && acchannels
		return out, sig_hwr, sig_hwr_ac
	elseif channels
		return out, sig_hwr
	elseif acchannels
		return out, sig_hwr_ac
	end
	return out
end



"""
	testsac(dur, f0, fb, fs; amps=[1.0 for k ∈ 1:12], 
           		    mistunings=[1.0 for k ∈ 1:12],
           		    phase=[0.0 for k ∈ 1:12],
           			actype=:fft,
           			compress = x -> x,
           			env = :none)

Generate the sac corresponding to a harmonic complex with mistuned harmonic(s).
"""
function testsac(dur, f0, fb, fs; amps=[1.0 for k ∈ 1:12], 
						    mistunings=[1.0 for k ∈ 1:12],
						    phase=[0.0 for k ∈ 1:12],
							actype=:fft,
							compress = x -> x,
							env = :none)
	sig = hc(dur + 0.05, f0, fs; amps=amps, mistunings=mistunings, phase=phase)
	y = sac(sig, fb, fs; actype=actype, dur=dur)
	y ./= maximum(y)
end

"""
	testsacfm(dur, f0, fb, fs, ffm, fp; amps=[1.0 for k ∈ 1:12], 
           		    mistunings=[1.0 for k ∈ 1:12],
           		    phase=[0.0 for k ∈ 1:12],
           			actype=:fft,
           			comp = x -> x,
           			env = :none)

Generate the sac corresponding to a harmonic complex with mistuned harmonic(s) with each component frequency modulated.
"""
function testsacfm(dur, f0, fb, fs, ffm, fp; amps=[1.0 for k ∈ 1:12], 
						    mistunings=[1.0 for k ∈ 1:12],
						    phase=[0.0 for k ∈ 1:12],
							actype=:fft,
							compress = x -> x,
							env = :none)
	sig = hcfm(dur + 0.05, f0, fs, ffm, fp; amps=amps, mistunings=mistunings, phase=phase)
	y = sac(sig, fb, fs; actype=actype, dur=dur)
	y ./= maximum(y)
end

"""
	gen_t(dur, fs)

Generate a vector of times for given duration and sampling rate.
"""
function gen_t(dur, fs)
	len = Int(round(dur*fs))
	collect(0:(len-1)) ./ fs
end


"""
	findpower2(n)
Find power of two greater than n.  (Used by edgestim for efficient FFT.)
"""
function findpower2(n)
	p = 1
	while 2 ^ p < n
		p += 1
	end
	return 2^p
end

"""
	edgestim(fe, fs, dur; kind=:hp, ftaper=16000, norm=true)

Generate an edge pitch stimulus in the frequency domain.  
"""
function edgestim(fe, fs, dur; kind=:hp, ftaper=16000, norm=true, seed=:none)
	if typeof(seed) <: Int 
		Random.seed!(seed)
	end
	len0 = Int(fs*dur)
	len = findpower2(len0)
	df = fs/len
	f = collect(range(0, step=df, length=len))
	mag = zeros(len)
	fe_idx = Int(round(fe/df))
	ft_idx = Int(round(ftaper/df))
	if kind == :hp
		mag[fe_idx:end-fe_idx] .= 1.0
		mid = Int(len/2)
		taper = collect(range(1.0, stop=10^(-2), length=mid-ft_idx))
		mag[ft_idx:(mid-1)] = taper
		mag[mid:(end-ft_idx-1)] = reverse(taper)
	elseif kind == :lp
		mag[1:fe_idx] .= 1.0
		mag[end-fe_idx:end] .= 1.0
	else
		error("Unknown kind of edge pitch. Select :hp or :lp")
	end
	phase = 2π .* rand(len)
	spec = mag .* exp.(im .* phase)
	sig = real.(ifft(spec))
	if norm
		e = energy(sig)
		sig ./ sqrt(e)
	end
	return sig[1:len0]
end

"""
	nepstim(fe, fs, dur_on, dur_off; kind=:hp, rampdur=0.005)

Edge pitch stimulus
"""
function nepstim(fe, fs, dur_on, dur_off; kind=:hp, rampdur=0.005, seed=:none)
	nep = edgestim(fe, fs, dur_on; kind=kind, seed=seed)	
	silence = zeros( Int(round(fs*dur_off)) )
	ramp = sinramp(rampdur, fs)
	nep[1:length(ramp)] .*= ramp
	nep[end-length(ramp)+1:end] .*= ramp[end:-1:1]
	[ silence; nep; silence ]
end

"""
	sinramp(dur, fs; kind=:on)

Sine ramp.
"""
function sinramp(dur, fs; kind=:on)
	t = collect( range(0, dur, step=1/fs) )
	f0 = 1.0 / (dur * 4.0)
	ramp = sin.( 2π*f0 .* t )
	if kind != :on
		ramp = ramp[end:-1:1]
	end
	return ramp
end



"""
	wn_lp(dur, fe, fs; ord=12, filttype=:lp)

Generate an edge pitch stimulus by filtering (less steep than edgestim).
"""
function wn_lp(dur, fe, fs; ord=12, filttype=:lp)
	len = Int(round(dur*fs))
	wn = randn(len)
	if filttype == :lp
		t = Lowpass(fe; fs=fs)
	else
		t = Highpass(fe; fs=fs)
	end
	m = Butterworth(ord)
	dfilt = digitalfilter(t, m)
	filt(dfilt, wn)
end

end
