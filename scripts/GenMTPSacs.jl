module GenMTPSacs

using Plots; pyplot()
using AuditoryFilters
using JLD

include("../util/EdgePitchUtils.jl")
include("../util/PeakPick.jl")
include("../util/StimGen.jl")
include("../util/DahlPlot.jl")
include("../core/Sieves.jl")
include("../core/Sigmas.jl")
include("../core/Salience.jl")
include("../core/Envelopes.jl")

# Trial Parameters (signals and kernel)
## Signal
fs = 48000
f0 = 155.0
mts = [0.887:0.0025:1.113;]
dur = 0.1
amps = [1.0 for k ∈ 1:12]
phase = [0.0 for k ∈ 1:12]
mtps = collect(1:5)

## Integration window
maxpeaks = Int(floor(Envelopes.hcc_maxdelay(1/f0) * f0))
maxpeaks = 4

## Auditory Filter Bank (if using full sacs)
fb = make_erb_filterbank(fs, 50, 50.0)

# Generate Stimuli
println("Generating Stimuli...")
sacs_p = []
@time for mtp ∈ mtps
	global sacs_p
	sacs = []
	ptunings = [1.0 for k ∈ 1:12]
	for mt ∈ mts
		ptunings[mtp] = mt
		sac = StimGen.testsac(dur, f0, fb, fs; amps=amps,
							  				   mistunings=ptunings,
											   phase=phase)
		push!(sacs, sac)
	end
	push!(sacs_p, sacs)
end
τ = collect(1:length(sacs_p[1][1])) ./ fs

save("/home/dahlbom/research/AdaptiveSieves/stimuli/MTPSacs155.jld",
	 "fs", fs,
	 "f0", f0,
	 "mts", mts,
	 "mtps", mtps,
	 "τ", τ,
	 "sacs_p", sacs_p)

end
