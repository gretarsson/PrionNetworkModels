#!/usr/bin/env julia

using CSV
using DataFrames
using Distributions
using PrionNetworkModels
using Random
using Statistics

function get_arg(flag::String, default::Union{Nothing,String}=nothing)
    idx = findfirst(==(flag), ARGS)
    isnothing(idx) && return default
    idx == length(ARGS) && error("Missing value for $flag")
    return ARGS[idx + 1]
end

function ks_statistic(x::AbstractVector, y::AbstractVector)
    xs = sort(collect(skipmissing(x)))
    ys = sort(collect(skipmissing(y)))
    n = length(xs)
    m = length(ys)
    (n == 0 || m == 0) && return NaN
    values = sort(unique(vcat(xs, ys)))
    i = 0
    j = 0
    d = 0.0
    for v in values
        while i < n && xs[i + 1] <= v
            i += 1
        end
        while j < m && ys[j + 1] <= v
            j += 1
        end
        d = max(d, abs(i / n - j / m))
    end
    return d
end

function ks_pvalue_asymptotic(d::Real, n::Integer, m::Integer)
    if !isfinite(d)
        return NaN
    end
    ne = n * m / (n + m)
    λ = (sqrt(ne) + 0.12 + 0.11 / sqrt(ne)) * d
    p = 2sum((-1)^(k - 1) * exp(-2k^2 * λ^2) for k in 1:100)
    return clamp(p, 0.0, 1.0)
end

function posterior_block(run_dir::AbstractString)
    posterior = load_posterior_hdf5(joinpath(run_dir, "posterior.h5"))
    spec = resolve_bundle_spec_paths(load_run_spec(joinpath(run_dir, "spec.toml")), run_dir)
    transport = build_transport_operator(spec.data.network; transport = spec.model.transport)
    priors = PrionNetworkModels.resolve_priors(spec, length(transport.labels))
    return posterior, spec, transport, priors
end

function param_draws(posterior, name::AbstractString)
    idx = findfirst(==(name), posterior.parameter_names)
    isnothing(idx) && return nothing
    return posterior.samples[:, idx]
end

function prior_draws(prior, n::Integer; seed::Integer=1)
    rng = MersenneTwister(seed)
    return rand(rng, prior, n)
end

function parameter_table(posterior, transport, priors, family::String; update_alpha::Float64=0.001)
    rows = DataFrame(
        region = String[],
        mean_post = Float64[],
        sd_post = Float64[],
        mean_prior = Float64[],
        sd_prior = Float64[],
        ks_pvalue = Float64[],
        updated = Int[],
    )

    labels = transport.labels
    n_draws = size(posterior.samples, 1)
    for (i, label) in enumerate(labels)
        name = "$(family)[$i]"
        draws = param_draws(posterior, name)
        isnothing(draws) && continue
        prior = getproperty(priors, Symbol(family))
        pdraws = prior_draws(prior, n_draws; seed = 10_000 + i)
        d = ks_statistic(draws, pdraws)
        p = ks_pvalue_asymptotic(d, length(draws), length(pdraws))
        push!(rows, (
            region = String(label),
            mean_post = mean(draws),
            sd_post = std(draws),
            mean_prior = mean(pdraws),
            sd_prior = std(pdraws),
            ks_pvalue = p,
            updated = isfinite(p) && p < update_alpha ? 1 : 0,
        ))
    end
    return rows
end

function main()
    run_dir = get_arg("--run", nothing)
    out_dir = get_arg("--out-dir", nothing)
    alpha = parse(Float64, get_arg("--update-alpha", "0.001"))
    isnothing(run_dir) && error("Usage: export_parameter_tables.jl --run runs/RUN_ID --out-dir paper-rf/results/parameters/RUN_ID")
    isnothing(out_dir) && error("Missing --out-dir")

    mkpath(out_dir)
    posterior, spec, transport, priors = posterior_block(run_dir)
    CSV.write(joinpath(out_dir, "beta.csv"), parameter_table(posterior, transport, priors, "beta"; update_alpha = alpha))
    CSV.write(joinpath(out_dir, "gamma.csv"), parameter_table(posterior, transport, priors, "gamma"; update_alpha = alpha))

    metadata = DataFrame(
        key = ["run_dir", "model", "transport", "n_draws", "n_regions", "update_alpha"],
        value = [abspath(run_dir), spec.model.name, spec.model.transport, string(size(posterior.samples, 1)), string(length(transport.labels)), string(alpha)],
    )
    CSV.write(joinpath(out_dir, "metadata.csv"), metadata)
    println(out_dir)
end

main()
