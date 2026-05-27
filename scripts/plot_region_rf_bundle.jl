#!/usr/bin/env julia

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using LinearAlgebra
using Plots
using PrionNetworkModels
using Statistics

function get_arg(flag::String, default::Union{Nothing,String}=nothing)
    idx = findfirst(==(flag), ARGS)
    isnothing(idx) && return default
    idx == length(ARGS) && error("Missing value for $flag")
    return ARGS[idx + 1]
end

function safe_slug(value::AbstractString)
    return replace(String(value), r"[^A-Za-z0-9_.-]+" => "_")
end

function diagnostics_paths(root::AbstractString)
    search_root = isdir(joinpath(root, "regional_runs")) ? joinpath(root, "regional_runs") : root
    paths = String[]
    for (dirpath, _, filenames) in walkdir(search_root)
        "diagnostics.csv" in filenames && push!(paths, joinpath(dirpath, "diagnostics.csv"))
    end
    sort!(paths)
    return paths
end

function posterior_summary_paths(root::AbstractString)
    search_root = isdir(joinpath(root, "regional_runs")) ? joinpath(root, "regional_runs") : root
    paths = String[]
    for (dirpath, _, filenames) in walkdir(search_root)
        if "diagnostics.csv" in filenames && "posterior_summary.csv" in filenames
            push!(paths, joinpath(dirpath, "posterior_summary.csv"))
        end
    end
    sort!(paths)
    return paths
end

function load_region_diagnostics(root::AbstractString)
    files = diagnostics_paths(root)
    isempty(files) && error("No diagnostics.csv files found under $root")
    tables = CSV.read.(files, Ref(DataFrame))
    df = reduce((a, b) -> vcat(a, b; cols = :union), tables)
    sort!(df, [:region_index])
    return df
end

function load_posterior_summaries(root::AbstractString)
    rows = DataFrame()
    for path in posterior_summary_paths(root)
        diag = CSV.read(joinpath(dirname(path), "diagnostics.csv"), DataFrame)
        posterior = CSV.read(path, DataFrame)
        posterior[!, :run_id] = fill(diag.run_id[1], nrow(posterior))
        posterior[!, :region_index] = fill(diag.region_index[1], nrow(posterior))
        posterior[!, :region] = fill(diag.region[1], nrow(posterior))
        posterior[!, :rank] = fill(diag.rank[1], nrow(posterior))
        append!(rows, posterior; cols = :union)
    end
    nrow(rows) == 0 && error("No posterior_summary.csv files found under $root")
    sort!(rows, [:region_index, :parameter])
    return rows
end

function region_rf_rhs!(du, u, p, t)
    alpha, beta, gamma = p
    x = u[1]
    y = u[2]
    du[1] = alpha * x * (beta - y - x)
    du[2] = gamma * x
    return nothing
end

function simulate_region_curve(timepoints, alpha, beta, gamma, u0; abstol = 1e-8, reltol = 1e-8)
    prob = ODEProblem(region_rf_rhs!, [u0, 0.0], (0.0, maximum(timepoints)), [alpha, beta, gamma])
    sol = solve(prob, Tsit5(); saveat = timepoints, abstol = abstol, reltol = reltol, maxiters = 50_000)
    return Array(sol)[1, :]
end

function observed_peak_ranking(obs)
    peaks = fill(-Inf, length(obs.labels))
    for i in eachindex(obs.labels)
        values = obs.mean[i, :]
        finite_values = values[isfinite.(values)]
        peaks[i] = isempty(finite_values) ? -Inf : maximum(finite_values)
    end
    ranked = sortperm(peaks; rev = true)
    ranks = similar(ranked)
    for (rank, region_idx) in enumerate(ranked)
        ranks[region_idx] = rank
    end
    finite_ranked = ranked[isfinite.(peaks[ranked])]
    return finite_ranked, ranks, peaks
end

function apply_observed_ranks!(diagnostics::DataFrame, obs)
    _, ranks, peaks = observed_peak_ranking(obs)
    diagnostics.observed_peak_mean = peaks[diagnostics.region_index]
    diagnostics.rank = ranks[diagnostics.region_index]
    sort!(diagnostics, [:region_index])
    return diagnostics
end

function assemble_predictions(diagnostics::DataFrame, obs)
    n_regions = length(obs.labels)
    n_times = length(obs.timepoints)
    pred = fill(NaN, n_regions, n_times)

    for row in eachrow(diagnostics)
        idx = Int(row.region_index)
        pred[idx, :] .= simulate_region_curve(
            obs.timepoints,
            row.alpha,
            row.beta,
            row.gamma,
            row.u0,
        )
    end

    return pred
end

function write_prediction_table(path::AbstractString, labels, timepoints, pred)
    df = DataFrame(region = labels)
    for (j, t) in enumerate(timepoints)
        df[!, string(t)] = pred[:, j]
    end
    CSV.write(path, df)
    return path
end

function r2_origin(observed, predicted)
    denom = sum(abs2, observed)
    denom == 0 && return NaN
    return 1 - sum(abs2, predicted .- observed) / denom
end

function predicted_observed_plot(obs, pred, output_path::AbstractString)
    x = vec(obs.mean)
    y = vec(pred)
    finite = isfinite.(x) .& isfinite.(y)
    x = x[finite]
    y = y[finite]
    minxy = min(minimum(x), minimum(y))
    maxxy = max(maximum(x), maximum(y))
    r2 = r2_origin(x, y)

    plt = scatter(
        x,
        y;
        xlabel = "Observed",
        ylabel = "Predicted",
        title = "Predicted vs Observed (R²=$(round(r2; digits = 3)))",
        legend = false,
        alpha = 0.65,
        markersize = 4,
        color = RGB(0 / 255, 71 / 255, 171 / 255),
        markerstrokecolor = :white,
        markerstrokewidth = 0.5,
        size = (720, 620),
    )
    plot!(plt, [minxy, maxxy], [minxy, maxxy]; color = :black, linestyle = :dash, linewidth = 2)
    savefig(plt, output_path)
    return output_path
end

function plot_rhat_diagnostics(posterior::DataFrame, output_dir::AbstractString)
    mkpath(output_dir)
    params = ["alpha", "beta", "gamma", "u0", "sigma"]
    filtered = posterior[in.(posterior.parameter, Ref(params)), :]
    CSV.write(joinpath(output_dir, "posterior_summary_long.csv"), posterior)
    CSV.write(joinpath(output_dir, "flagged_parameters.csv"), filtered[coalesce.(filtered.rhat .> 1.05, false), :])

    subplots = Plots.Plot[]
    for param in params
        sub = filtered[filtered.parameter .== param, :]
        plt = scatter(
            sub.rank,
            sub.rhat;
            xlabel = "Observed pathology rank",
            ylabel = "R̂",
            title = param,
            legend = false,
            markersize = 3,
            alpha = 0.75,
            color = RGB(0 / 255, 71 / 255, 171 / 255),
        )
        hline!(plt, [1.01, 1.05, 1.10]; color = [:gray :orange :red], linestyle = :dash, linewidth = 1.4)
        push!(subplots, plt)
    end
    panel = plot(subplots...; layout = (3, 2), size = (1100, 1050), plot_title = "Region-wise RF R̂ by Parameter")
    savefig(panel, joinpath(output_dir, "rhat_by_rank.pdf"))
    savefig(panel, joinpath(output_dir, "rhat_by_rank.png"))

    histplots = Plots.Plot[]
    for param in params
        sub = filtered[filtered.parameter .== param, :]
        plt = histogram(
            sub.rhat;
            xlabel = "R̂",
            ylabel = "Regions",
            title = param,
            legend = false,
            bins = 40,
            color = RGB(0 / 255, 71 / 255, 171 / 255),
            alpha = 0.75,
        )
        vline!(plt, [1.01, 1.05, 1.10]; color = [:gray :orange :red], linestyle = :dash, linewidth = 1.4)
        push!(histplots, plt)
    end
    hist_panel = plot(histplots...; layout = (3, 2), size = (1100, 1050), plot_title = "Region-wise RF R̂ Distributions")
    savefig(hist_panel, joinpath(output_dir, "rhat_histograms.pdf"))
    savefig(hist_panel, joinpath(output_dir, "rhat_histograms.png"))

    rhat_summary = combine(groupby(filtered, :parameter),
        :rhat => length => :n_regions,
        :rhat => maximum => :max_rhat,
        :rhat => (x -> mean(x .> 1.01)) => :frac_gt_1_01,
        :rhat => (x -> mean(x .> 1.05)) => :frac_gt_1_05,
        :rhat => (x -> mean(x .> 1.10)) => :frac_gt_1_10,
    )
    CSV.write(joinpath(output_dir, "rhat_summary.csv"), rhat_summary)
    return output_dir
end

function plot_parameter_maps(diagnostics::DataFrame, output_dir::AbstractString)
    mkpath(output_dir)
    params = ["alpha", "beta", "gamma", "u0", "sigma"]
    subplots = Plots.Plot[]
    for param in params
        plt = scatter(
            diagnostics.rank,
            diagnostics[!, param];
            xlabel = "Observed pathology rank",
            ylabel = param,
            title = param,
            legend = false,
            markersize = 3,
            alpha = 0.75,
            color = RGB(196 / 255, 54 / 255, 22 / 255),
        )
        push!(subplots, plt)
    end
    panel = plot(subplots...; layout = (3, 2), size = (1100, 1050), plot_title = "Posterior Means by Observed Pathology Rank")
    savefig(panel, joinpath(output_dir, "posterior_means_by_rank.pdf"))
    savefig(panel, joinpath(output_dir, "posterior_means_by_rank.png"))
    return output_dir
end

function smooth_region_predictions(diagnostics::DataFrame, obs; n_dense::Int = 300)
    dense_time = collect(range(0.0, maximum(obs.timepoints); length = n_dense))
    n_regions = length(obs.labels)
    mean_path = fill(NaN, n_regions, n_dense)
    lower50 = similar(mean_path)
    upper50 = similar(mean_path)
    lower90 = similar(mean_path)
    upper90 = similar(mean_path)

    inner_z = quantile(Normal(), 0.75)
    outer_z = quantile(Normal(), 0.95)

    for row in eachrow(diagnostics)
        idx = Int(row.region_index)
        pred = simulate_region_curve(dense_time, row.alpha, row.beta, row.gamma, row.u0)
        sigma = row.sigma
        mean_path[idx, :] .= pred
        lower50[idx, :] .= max.(pred .- inner_z * sigma, 0.0)
        upper50[idx, :] .= pred .+ inner_z * sigma
        lower90[idx, :] .= max.(pred .- outer_z * sigma, 0.0)
        upper90[idx, :] .= pred .+ outer_z * sigma
    end

    return (timepoints = dense_time, mean = mean_path, lower50 = lower50, upper50 = upper50, lower90 = lower90, upper90 = upper90)
end

function retrodiction_plot(output_path::AbstractString, obs, smooth, region_idx::Integer; ymax::Float64)
    observed_color = RGB(0 / 255, 71 / 255, 171 / 255)
    plt = plot(
        xlabel = "Time",
        ylabel = "Pathology",
        title = obs.labels[region_idx],
        legend = :topleft,
        ylims = (0.0, ymax),
        linewidth = 3,
        size = (720, 420),
    )
    plot!(
        plt,
        smooth.timepoints,
        smooth.mean[region_idx, :];
        ribbon = (
            smooth.mean[region_idx, :] .- smooth.lower90[region_idx, :],
            smooth.upper90[region_idx, :] .- smooth.mean[region_idx, :],
        ),
        fillalpha = 0.18,
        fillcolor = :gray70,
        linealpha = 0.0,
        label = "90% noise band",
    )
    plot!(
        plt,
        smooth.timepoints,
        smooth.mean[region_idx, :];
        ribbon = (
            smooth.mean[region_idx, :] .- smooth.lower50[region_idx, :],
            smooth.upper50[region_idx, :] .- smooth.mean[region_idx, :],
        ),
        fillalpha = 0.28,
        fillcolor = :gray45,
        linealpha = 0.0,
        label = "50% noise band",
    )
    plot!(plt, smooth.timepoints, smooth.mean[region_idx, :]; color = :black, linewidth = 3, label = "Posterior mean")
    scatter!(
        plt,
        obs.timepoints,
        obs.mean[region_idx, :];
        yerror = obs.se[region_idx, :],
        label = "Observed mean ± SE",
        color = observed_color,
        markersize = 5,
        markerstrokecolor = :white,
        markerstrokewidth = 0.7,
    )
    savefig(plt, output_path)
    return output_path
end

function retrodiction_panels(obs, smooth, output_dir::AbstractString; n_panels::Int = 3, regions_per_panel::Int = 4)
    mkpath(output_dir)
    ranked, _, peaks = observed_peak_ranking(obs)
    total_regions = min(length(ranked), n_panels * regions_per_panel)
    observed_color = RGB(0 / 255, 71 / 255, 171 / 255)

    for panel_idx in 1:ceil(Int, total_regions / regions_per_panel)
        start_idx = (panel_idx - 1) * regions_per_panel + 1
        stop_idx = min(panel_idx * regions_per_panel, total_regions)
        region_indices = ranked[start_idx:stop_idx]
        ymax = maximum(skipmissing(vec(obs.mean[region_indices, :])))
        ymax = max(ymax, maximum(smooth.upper90[region_indices, :]))
        ymax = max(ymax, 1e-3) * 1.08

        subplots = Plots.Plot[]
        for (rank_idx, region_idx) in enumerate(region_indices)
            rank = start_idx + rank_idx - 1
            plt = plot(
                xlabel = "Time",
                ylabel = "Pathology",
                title = "#$rank $(obs.labels[region_idx])",
                legend = rank_idx == 1 ? :topleft : false,
                ylims = (0.0, ymax),
                linewidth = 2.5,
            )
            plot!(
                plt,
                smooth.timepoints,
                smooth.mean[region_idx, :];
                ribbon = (
                    smooth.mean[region_idx, :] .- smooth.lower90[region_idx, :],
                    smooth.upper90[region_idx, :] .- smooth.mean[region_idx, :],
                ),
                fillalpha = 0.16,
                fillcolor = :gray70,
                linealpha = 0.0,
                label = "90% noise band",
            )
            plot!(plt, smooth.timepoints, smooth.mean[region_idx, :]; color = :black, linewidth = 2.7, label = "Posterior mean")
            scatter!(
                plt,
                obs.timepoints,
                obs.mean[region_idx, :];
                yerror = obs.se[region_idx, :],
                label = "Observed mean ± SE",
                color = observed_color,
                markersize = 4,
                markerstrokecolor = :white,
                markerstrokewidth = 0.6,
            )
            push!(subplots, plt)
        end

        panel = plot(
            subplots...;
            layout = length(region_indices) <= 2 ? (1, length(region_indices)) : (2, 2),
            size = (1100, 820),
            plot_title = "Top observed pathology regions $(start_idx)-$(stop_idx)",
        )
        filename = "top_observed_pathology_$(start_idx)_to_$(stop_idx)"
        savefig(panel, joinpath(output_dir, "$filename.pdf"))
        savefig(panel, joinpath(output_dir, "$filename.png"))
    end

    ranking = DataFrame(
        rank = collect(1:length(ranked)),
        region_index = ranked,
        region = obs.labels[ranked],
        peak_observed_mean = peaks[ranked],
    )
    CSV.write(joinpath(output_dir, "top_observed_pathology_region_ranking.csv"), ranking)
    return output_dir
end

function retrodiction_plots(obs, diagnostics::DataFrame, output_dir::AbstractString)
    mkpath(output_dir)
    smooth = smooth_region_predictions(diagnostics, obs)

    global_ymax = maximum(skipmissing(vec(obs.mean)))
    global_ymax = max(global_ymax, maximum(smooth.upper90))
    global_ymax = max(global_ymax, 1e-3) * 1.05

    for i in eachindex(obs.labels)
        retrodiction_plot(
            joinpath(output_dir, "retrodiction_$(safe_slug(obs.labels[i])).pdf"),
            obs,
            smooth,
            i;
            ymax = global_ymax,
        )
    end
    retrodiction_panels(obs, smooth, joinpath(output_dir, "top_pathology_panels"))
    return output_dir
end

function main()
    root = get_arg("--root", nothing)
    observations = get_arg("--observations", nothing)
    network = get_arg("--network", nothing)
    isnothing(root) && error("Usage: plot_region_rf_bundle.jl --root runs/region_rf/DATASET --observations observations.csv --network network.csv [--out plots_dir]")
    isnothing(observations) && error("Missing --observations")
    isnothing(network) && error("Missing --network")
    out_dir = get_arg("--out", joinpath(root, "plots"))

    diagnostics = load_region_diagnostics(root)
    posterior = load_posterior_summaries(root)
    obs = let
        pathology = process_pathology(observations; network_csv = network)
        summary = summarize_over_replicates(pathology.data)
        (
            labels = pathology.labels,
            timepoints = pathology.timepoints,
            mean = summary.mean,
            sd = summary.sd,
            se = summary.se,
            n = summary.n,
        )
    end
    apply_observed_ranks!(diagnostics, obs)
    rank_map = Dict(row.region_index => row.rank for row in eachrow(diagnostics))
    if "rank" in names(posterior)
        posterior.rank = [rank_map[row.region_index] for row in eachrow(posterior)]
    end

    mkpath(out_dir)
    CSV.write(joinpath(root, "region_rf_summary.csv"), diagnostics)
    CSV.write(joinpath(root, "region_rf_posterior_summary_long.csv"), posterior)

    pred = assemble_predictions(diagnostics, obs)
    write_prediction_table(joinpath(root, "predictions_train.csv"), obs.labels, obs.timepoints, pred)
    write_prediction_table(joinpath(out_dir, "predictions_train.csv"), obs.labels, obs.timepoints, pred)

    predicted_observed_plot(obs, pred, joinpath(out_dir, "predicted_vs_observed.pdf"))
    predicted_observed_plot(obs, pred, joinpath(out_dir, "predicted_vs_observed.png"))
    plot_rhat_diagnostics(posterior, joinpath(out_dir, "diagnostics"))
    plot_parameter_maps(diagnostics, joinpath(out_dir, "diagnostics"))
    retrodiction_plots(obs, diagnostics, joinpath(out_dir, "retrodiction"))

    println(out_dir)
end

main()
