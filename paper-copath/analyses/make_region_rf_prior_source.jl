#!/usr/bin/env julia

using CSV
using DataFrames
using HDF5
using Random
using Statistics

const PARAMETERS = ["alpha", "beta", "gamma"]

function get_arg(flag::String, default::Union{Nothing,String}=nothing)
    idx = findfirst(==(flag), ARGS)
    isnothing(idx) && return default
    idx == length(ARGS) && error("Missing value for $flag")
    return ARGS[idx + 1]
end

function prior_source_from_region_rf(adjusted_dir::AbstractString, dataset::AbstractString, out_path::AbstractString; n_draws::Int=1000, seed::Int=1, min_sd::Float64=1e-6)
    posterior_path = joinpath(adjusted_dir, "$(dataset)_region_rf_posterior_summary_long.csv")
    isfile(posterior_path) || error("Missing REGION-RF posterior summary: $posterior_path")

    posterior = CSV.read(posterior_path, DataFrame)
    posterior.parameter = String.(posterior.parameter)
    posterior = posterior[in.(posterior.parameter, Ref(PARAMETERS)), :]

    region_indices = sort(unique(Int.(posterior.region_index)))
    N = maximum(region_indices)
    region_indices == collect(1:N) || error("Expected contiguous region_index values 1:$N in $posterior_path")

    region_labels = Vector{String}(undef, N)
    parameter_names = String[]
    means = Float64[]
    sds = Float64[]
    for parameter in PARAMETERS
        for region_idx in 1:N
            rows = posterior[(posterior.parameter .== parameter) .& (posterior.region_index .== region_idx), :]
            nrow(rows) == 1 || error("Expected one row for $parameter[$region_idx], found $(nrow(rows))")
            if parameter == first(PARAMETERS)
                region_labels[region_idx] = String(rows.region[1])
            end
            push!(parameter_names, "$parameter[$region_idx]")
            push!(means, Float64(rows.mean[1]))
            sd_col = :std in propertynames(rows) ? :std : :sd
            sd = Float64(rows[1, sd_col])
            push!(sds, isfinite(sd) && sd > 0 ? sd : min_sd)
        end
    end

    rng = MersenneTwister(seed)
    samples = Matrix{Float64}(undef, n_draws, length(parameter_names))
    for j in axes(samples, 2)
        samples[:, j] .= means[j] .+ sds[j] .* randn(rng, n_draws)
    end

    mkpath(dirname(out_path))
    h5open(out_path, "w") do h5
        h5["chains/samples"] = samples
        h5["chains/parameter_names"] = parameter_names
        h5["data/region_labels"] = region_labels
        h5["spec/model_name"] = "REGION-RF"
        h5["spec/source_dataset"] = dataset
    end
    return out_path
end

function main()
    adjusted_dir = get_arg("--adjusted-region-rf-dir", "paper-copath/results/region_rf_iterative_drop_low_likelihood_chains")
    dataset = get_arg("--dataset", nothing)
    isnothing(dataset) && error("Usage: make_region_rf_prior_source.jl --dataset syn_mapt [--out path.h5]")
    out_path = get_arg("--out", joinpath("paper-copath/results/diff_rf_regional_priors", "$(dataset)_region_rf_prior_source.h5"))
    n_draws = parse(Int, get_arg("--n-draws", "1000"))
    seed = parse(Int, get_arg("--seed", "1"))
    min_sd = parse(Float64, get_arg("--min-sd", "1e-6"))
    println(prior_source_from_region_rf(adjusted_dir, dataset, out_path; n_draws=n_draws, seed=seed, min_sd=min_sd))
end

main()
