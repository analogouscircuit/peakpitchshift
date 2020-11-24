module ACUtils

using FFTW

export ac_true, ac_fft

function ac_true(sig; norm=true)
      out = zeros(length(sig))
      for k ∈ 1:length(out)
            out[k] = sig[k:end]' * sig[1:(end-k+1)]
      end
      if norm
            return out ./ maximum(out)
      end
      return out
end

function ac_fft(sig; norm=true)
	sig_padded = vcat(sig, zeros(length(sig)))
	spec = fft(sig_padded)
	out = real.(ifft( conj(spec) .* spec ))
	out = out[1:length(sig)]
	if norm
		return out ./ maximum(out)
	end
	return out
end

end
