module Sigmas

function uniform(;σ0=5e-5, num_bumps=6)
	return [σ0 for k in 1:num_bumps]
end

function uniform_scaled(τ0; σ0=5e-5, num_bumps=6)
	σ = σ0 * (τ0*155.0)	# relative to width at 155.0 Hz (Darwin and Ciocca)
	return [σ for k in 1:num_bumps]
end

function nonuniform_scaled(τ0; σ0=5e-5, num_bumps=6)
	σ = σ0 * (τ0*155.0)	# relative to width at 155.0 Hz (Darwin and Ciocca)
	return [σ*k for k in 1:num_bumps]
end

function σ0_template(τ0; κ=0.07168)
	 κ * τ0^2
end

function nonuniform_envinv(τ0, env; σ0=0.0, num_bumps=6) 
	envpts = env([τ0 for p in 1:num_bumps], τ0) ./ env([τ0*p for p in 1:num_bumps], τ0)
	sigmas = envpts .* σ0_template(τ0) 
	sigmas ./ [p*τ0 for p in 1:num_bumps]
end



end
