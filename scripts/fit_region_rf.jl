#!/usr/bin/env julia

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using HDF5
using LinearAlgebra
using MCMCChains
using Plots
using PrionNetworkModels
using Random
using ReverseDiff
using SciMLSensitivity
using Statistics
using Turing

function get_arg(flag::String, default::Union{Nothing,String}=nothing)
    idx = findfirst(==(flag), ARGS)
    isnothing(idx) && return default
    idx == length(ARGS) && error("Missing value for $flag")
    return ARGS[idx + 1]
end

has_flag(flag::String) = flag in ARGS

function parse_int_list(value::AbstractString)
    return parse.(Int, [strip(part) for part in split(value, ",") if !isempty(strip(part))])
end

function safe_slug(value::AbstractString)
    return replace(String(value), r"[^A-Za-z0-9_.-]+" => "_")
end

function observed_peak_ranking(pathology)
    summary = summarize_over_replicates(pathology.data)
    peaks = fill(-Inf, length(pathology.labels))
    for i in eachindex(pathology.labels)
        values = summary.mean[i, :]
        finite_values = values[isfinite.(values)]
        peaks[i] = isempty(finite_values) ? -Inf : maximum(finite_values)
    end
    ranked = sortperm(peaks; rev = true)
    ranks = similar(ranked)
    for (rank, region_idx) in enumerate(ranked)
        ranks[region_idx] = rank
    end
    return ranked, ranks, peaks, summary
end

function region_observations(data, region_idx::Integer; mean_data::Bool=false)
    region_data = data[region_idx:region_idx, :, :]
    obs = PrionNetworkModels.finite_observations(region_data; mean_data = mean_data)
    if obs.mean_data
        time_indices = findall(vec(obs.mask))
    else
        _, n_times, n_samples = size(region_data)
        time_indices = Int[]
        for linear_idx in obs.nonmissing
            _, t, _ = Tuple(CartesianIndices((1, n_times, n_samples))[linear_idx])
            push!(time_indices, t)
        end
    end
    return (values = obs.values, time_indices = time_indices, mean_data = obs.mean_data)
end

function region_rf_rhs!(du, u, p, t)
    alpha, beta, gamma = p
    x = u[1]
    y = u[2]
    du[1] = alpha * x * (beta - y - x)
    du[2] = gamma * x
    return nothing
end

@model function region_rf_model(timepoints::Vector{Float64}, obs_vals::Vector{Float64}, obs_time_indices::Vector{Int}, maxiters::Int, u0_prior_sd::Float64, alpha_prior_sd::Float64)
    alpha ~ truncated(Normal(0.0, alpha_prior_sd), 0.0, Inf)
    beta ~ Normal(0.0, 1.0)
    gamma ~ truncated(Normal(0.0, 0.1), 0.0, Inf)
    u0 ~ truncated(Normal(0.0, u0_prior_sd), 0.0, Inf)
    sigma ~ LogNormal(0, 1)

    prob = ODEProblem(region_rf_rhs!, [u0, 0.0], (0.0, maximum(timepoints)), [alpha, beta, gamma])
    sol = solve(
        prob,
        Tsit5();
        saveat = timepoints,
        sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP(true)),
        abstol = 1e-6,
        reltol = 1e-3,
        maxiters = maxiters,
    )
    pred = Array(sol)[1, :]
    pred_vec = pred[obs_time_indices]
    obs_vals ~ MvNormal(pred_vec, sigma^2 * I)
end

function posterior_summary(chain::Chains)
    df = DataFrame(MCMCChains.summarystats(chain))
    rename!(df, :parameters => :parameter)
    return df
end

function posterior_mean(chain::Chains, name::AbstractString)
    return mean(vec(Array(chain[Symbol(name)])))
end

function chain_array(chain::Chains)
    per_chain = Array(chain; append_chains = false)
    n_chain = length(per_chain)
    n_iter, n_param = size(first(per_chain))
    samples = Array{Float64}(undef, n_iter, n_param, n_chain)
    for chain_idx in 1:n_chain
        samples[:, :, chain_idx] = per_chain[chain_idx]
    end
    return samples
end

function write_chain_h5(path::AbstractString, chain::Chains)
    samples = chain_array(chain)
    param_names = String.(names(chain, :parameters))
    n_iter, n_param, n_chain = size(samples)
    flat = reshape(permutedims(samples, (1, 3, 2)), n_iter * n_chain, n_param)
    chain_ids = reduce(vcat, [fill(c, n_iter) for c in 1:n_chain])
    h5open(path, "w") do h5
        h5["chains/samples"] = flat
        h5["chains/parameter_names"] = param_names
        h5["chains/chain_ids"] = chain_ids
    end
end

function write_trace_plots(output_dir::AbstractString, chain::Chains)
    mkpath(output_dir)
    samples = chain_array(chain)
    param_names = String.(names(chain, :parameters))
    colors = [
        RGB(0 / 255, 71 / 255, 171 / 255),
        RGB(196 / 255, 54 / 255, 22 / 255),
        RGB(0 / 255, 136 / 255, 55 / 255),
        RGB(123 / 255, 31 / 255, 162 / 255),
    ]

    for (param_idx, param_name) in enumerate(param_names)
        plt = plot(
            xlabel = "Iteration",
            ylabel = param_name,
            title = "Trace: $param_name",
            legend = :topright,
            size = (860, 480),
        )
        for chain_idx in axes(samples, 3)
            plot!(
                plt,
                axes(samples, 1),
                samples[:, param_idx, chain_idx];
                label = "Chain $chain_idx",
                linewidth = 1.6,
                alpha = 0.9,
                color = colors[mod1(chain_idx, length(colors))],
            )
        end
        savefig(plt, joinpath(output_dir, "trace_$(safe_slug(param_name)).pdf"))
        savefig(plt, joinpath(output_dir, "trace_$(safe_slug(param_name)).png"))
    end
    return output_dir
end

function plot_region_fit(output_path::AbstractString, chain::Chains, timepoints, obs_mean, obs_se, region_label::AbstractString)
    params = [posterior_mean(chain, "alpha"), posterior_mean(chain, "beta"), posterior_mean(chain, "gamma")]
    u0 = posterior_mean(chain, "u0")
    sigma = posterior_mean(chain, "sigma")
    dense_timepoints = collect(range(0.0, maximum(timepoints); length = 300))
    prob = ODEProblem(region_rf_rhs!, [u0, 0.0], (0.0, maximum(timepoints)), params)
    sol = solve(prob, Tsit5(); saveat = dense_timepoints, abstol = 1e-8, reltol = 1e-8)
    pred = Array(sol)[1, :]

    upper90 = pred .+ quantile(Normal(), 0.95) * sigma
    lower90 = max.(pred .- quantile(Normal(), 0.95) * sigma, 0.0)

    plt = plot(
        dense_timepoints,
        pred;
        ribbon = (pred .- lower90, upper90 .- pred),
        fillalpha = 0.16,
        fillcolor = :gray70,
        color = :black,
        linewidth = 2.8,
        label = "Posterior mean",
        xlabel = "Time",
        ylabel = "Pathology",
        title = region_label,
        size = (760, 520),
    )
    scatter!(
        plt,
        timepoints,
        obs_mean;
        yerror = obs_se,
        label = "Observed mean ± SE",
        color = RGB(0 / 255, 71 / 255, 171 / 255),
        markersize = 5,
        markerstrokecolor = :white,
        markerstrokewidth = 0.7,
    )
    savefig(plt, output_path)
    return output_path
end

function fit_region(pathology, summary, region_idx::Integer; rank::Integer, output_root::AbstractString, run_prefix::AbstractString, n_samples::Int, n_warmup::Int, n_chains::Int, mean_data::Bool, progress::Bool, maxiters::Int, seed::Int, u0_prior_sd::Float64, alpha_prior_sd::Float64, write_traces::Bool)
    label = pathology.labels[region_idx]
    run_id = "$(run_prefix)_rank$(rank)_region$(region_idx)_$(safe_slug(label))"
    run_dir = joinpath(output_root, run_id)
    mkpath(run_dir)

    obs = region_observations(pathology.data, region_idx; mean_data = mean_data)
    model = region_rf_model(Float64.(pathology.timepoints), obs.values, obs.time_indices, maxiters, u0_prior_sd, alpha_prior_sd)

    start_time = time()
    chains = Chains[]
    for chain_idx in 1:n_chains
        Random.seed!(seed + 10_000 * region_idx + chain_idx)
        println("  Chain $chain_idx/$n_chains with seed $(seed + 10_000 * region_idx + chain_idx)")
        push!(
            chains,
            sample(
                model,
                NUTS(n_warmup, 0.65; adtype = AutoReverseDiff()),
                n_samples;
                progress = progress,
            ),
        )
    end
    chain = length(chains) == 1 ? first(chains) : chainscat(first(chains), chains[2:end]...)
    elapsed = time() - start_time
    samples = chain_array(chain)
    saved_iterations, _, saved_chains = size(samples)

    CSV.write(joinpath(run_dir, "posterior_summary.csv"), posterior_summary(chain))
    write_chain_h5(joinpath(run_dir, "posterior.h5"), chain)
    write_traces && write_trace_plots(joinpath(run_dir, "trace"), chain)

    pred_df = DataFrame(region = [label])
    params = [posterior_mean(chain, "alpha"), posterior_mean(chain, "beta"), posterior_mean(chain, "gamma")]
    u0 = posterior_mean(chain, "u0")
    prob = ODEProblem(region_rf_rhs!, [u0, 0.0], (0.0, maximum(pathology.timepoints)), params)
    sol = solve(prob, Tsit5(); saveat = pathology.timepoints, abstol = 1e-8, reltol = 1e-8)
    pred = Array(sol)[1, :]
    for (j, t) in enumerate(pathology.timepoints)
        pred_df[!, string(t)] = [pred[j]]
    end
    CSV.write(joinpath(run_dir, "predictions_train.csv"), pred_df)

    plot_region_fit(
        joinpath(run_dir, "fit.pdf"),
        chain,
        pathology.timepoints,
        summary.mean[region_idx, :],
        summary.se[region_idx, :],
        label,
    )
    plot_region_fit(
        joinpath(run_dir, "fit.png"),
        chain,
        pathology.timepoints,
        summary.mean[region_idx, :],
        summary.se[region_idx, :],
        label,
    )

    diagnostics = DataFrame(
        run_id = [run_id],
        region_index = [region_idx],
        region = [label],
        rank = [rank],
        elapsed_seconds = [elapsed],
        n_samples = [n_samples],
        n_warmup = [n_warmup],
        n_chains = [n_chains],
        maxiters = [maxiters],
        seed = [seed],
        u0_prior_sd = [u0_prior_sd],
        alpha_prior_sd = [alpha_prior_sd],
        saved_iterations_per_chain = [saved_iterations],
        saved_chains = [saved_chains],
        saved_draws_total = [saved_iterations * saved_chains],
        alpha = [posterior_mean(chain, "alpha")],
        beta = [posterior_mean(chain, "beta")],
        gamma = [posterior_mean(chain, "gamma")],
        u0 = [posterior_mean(chain, "u0")],
        sigma = [posterior_mean(chain, "sigma")],
    )
    CSV.write(joinpath(run_dir, "diagnostics.csv"), diagnostics)
    println("Finished $run_id in $(round(elapsed; digits=1)) seconds")
    return diagnostics
end

function main()
    root = dirname(@__DIR__)
    observations = get_arg("--observations", joinpath(root, "paper-rf", "data", "striatum", "observations.csv"))
    network = get_arg("--network", joinpath(root, "paper-rf", "data", "striatum", "network.csv"))
    output_root = get_arg("--out-root", joinpath(root, "runs", "region_rf"))
    run_prefix = get_arg("--run-prefix", "striatum_REGION-RF")
    top_n = parse(Int, get_arg("--top", "4"))
    regions_arg = get_arg("--regions", nothing)
    region_index_arg = get_arg("--region-index", nothing)
    n_samples = parse(Int, get_arg("--samples", "1000"))
    n_warmup = parse(Int, get_arg("--warmup", "1000"))
    n_chains = parse(Int, get_arg("--chains", "1"))
    maxiters = parse(Int, get_arg("--maxiters", "10000"))
    seed = parse(Int, get_arg("--seed", "8675309"))
    u0_prior_sd = parse(Float64, get_arg("--u0-prior-sd", "0.01"))
    alpha_prior_sd = parse(Float64, get_arg("--alpha-prior-sd", "0.1"))
    mean_data = has_flag("--mean-data")
    progress = has_flag("--progress")
    no_root_summary = has_flag("--no-root-summary")
    skip_traces = has_flag("--skip-traces")

    pathology = process_pathology(observations; network_csv = network)
    ranked, ranks, peaks, summary = observed_peak_ranking(pathology)
    selected = if !isnothing(region_index_arg)
        [parse(Int, region_index_arg)]
    elseif !isnothing(regions_arg)
        parse_int_list(regions_arg)
    else
        ranked[1:min(top_n, length(ranked))]
    end
    mkpath(output_root)

    rows = DataFrame()
    for (position, region_idx) in enumerate(selected)
        rank = ranks[region_idx]
        println("Fitting rank $rank region $region_idx ($(pathology.labels[region_idx])); peak observed mean = $(peaks[region_idx])")
        diagnostics = fit_region(
            pathology,
            summary,
            region_idx;
            rank = rank,
            output_root = output_root,
            run_prefix = run_prefix,
            n_samples = n_samples,
            n_warmup = n_warmup,
            n_chains = n_chains,
            mean_data = mean_data,
            progress = progress,
            maxiters = maxiters,
            seed = seed,
            u0_prior_sd = u0_prior_sd,
            alpha_prior_sd = alpha_prior_sd,
            write_traces = !skip_traces,
        )
        append!(rows, diagnostics)
        no_root_summary || CSV.write(joinpath(output_root, "$(run_prefix)_summary.csv"), rows)
    end

    no_root_summary || println(joinpath(output_root, "$(run_prefix)_summary.csv"))
end

main()
