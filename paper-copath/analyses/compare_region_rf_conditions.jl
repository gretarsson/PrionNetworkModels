#!/usr/bin/env julia

using CSV
using DataFrames
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

function ensure_summary(root::AbstractString)
    path = joinpath(root, "region_rf_summary.csv")
    isfile(path) || error("Missing REGION-RF summary: $path. Run paper-copath/plot_region_rf_copath.sh first.")
    return CSV.read(path, DataFrame)
end

function ensure_posterior_summary(root::AbstractString)
    path = joinpath(root, "region_rf_posterior_summary_long.csv")
    isfile(path) || error("Missing REGION-RF posterior summary: $path. Run paper-copath/plot_region_rf_copath.sh first.")
    return CSV.read(path, DataFrame)
end

function pathology_peaks(observations::AbstractString, network::AbstractString)
    pathology = process_pathology(observations; network_csv = network)
    summary = summarize_over_replicates(pathology.data)
    peaks = fill(NaN, length(pathology.labels))
    for i in eachindex(pathology.labels)
        values = summary.mean[i, :]
        finite_values = values[isfinite.(values)]
        peaks[i] = isempty(finite_values) ? NaN : maximum(finite_values)
    end
    return DataFrame(region_index = collect(eachindex(pathology.labels)), region = pathology.labels, peak = peaks)
end

function finite_pair_max(a, b)
    vals = Float64[]
    isfinite(a) && push!(vals, Float64(a))
    isfinite(b) && push!(vals, Float64(b))
    return isempty(vals) ? NaN : maximum(vals)
end

function parameter_rhats(root::AbstractString, suffix::AbstractString)
    posterior = ensure_posterior_summary(root)
    keep = in.(posterior.parameter, Ref(["alpha", "beta", "gamma"]))
    posterior = posterior[keep, [:region_index, :parameter, :rhat]]
    rename!(posterior, :rhat => Symbol("rhat_$suffix"))
    wide = unstack(posterior, :region_index, :parameter, Symbol("rhat_$suffix"))
    rename!(wide, Dict(name => Symbol("$(name)_rhat_$suffix") for name in names(wide) if name != "region_index"))
    return wide
end

function comparison_table(project_root::AbstractString, protein::AbstractString)
    app_root = joinpath(project_root, "runs", "region_rf", "copath_$(protein)_app")
    mapt_root = joinpath(project_root, "runs", "region_rf", "copath_$(protein)_mapt")
    network = joinpath(project_root, "paper-copath", "data", "network.csv")
    app_obs = joinpath(project_root, "paper-copath", "data", "$(protein)_pathology_app.csv")
    mapt_obs = joinpath(project_root, "paper-copath", "data", "$(protein)_pathology_mapt.csv")

    app = ensure_summary(app_root)
    mapt = ensure_summary(mapt_root)
    app = app[:, [:region_index, :region, :rank, :alpha, :beta, :gamma, :u0, :sigma]]
    mapt = mapt[:, [:region_index, :region, :rank, :alpha, :beta, :gamma, :u0, :sigma]]
    rename!(app, Dict(name => Symbol("$(name)_app") for name in names(app) if !(name in ["region_index", "region"])))
    rename!(mapt, Dict(name => Symbol("$(name)_mapt") for name in names(mapt) if !(name in ["region_index", "region"])))
    df = innerjoin(app, mapt, on = [:region_index, :region])

    peaks_app = rename(pathology_peaks(app_obs, network), :peak => :peak_app)
    peaks_mapt = rename(pathology_peaks(mapt_obs, network), :peak => :peak_mapt)
    df = leftjoin(df, peaks_app, on = [:region_index, :region])
    df = leftjoin(df, peaks_mapt, on = [:region_index, :region])
    df.peak_any = finite_pair_max.(df.peak_app, df.peak_mapt)
    df.active_any = isfinite.(df.peak_any) .& (df.peak_any .> 0)
    rank_values = ifelse.(isfinite.(df.peak_any), df.peak_any, -Inf)
    rank_idx = sortperm(rank_values; rev = true)
    condition_rank = similar(df.region_index)
    for (rank, idx) in enumerate(rank_idx)
        condition_rank[idx] = rank
    end
    df.condition_pair_rank = condition_rank

    df = leftjoin(df, parameter_rhats(app_root, "app"), on = :region_index)
    df = leftjoin(df, parameter_rhats(mapt_root, "mapt"), on = :region_index)

    for param in ["alpha", "beta", "gamma"]
        app_col = Symbol("$(param)_app")
        mapt_col = Symbol("$(param)_mapt")
        diff_col = Symbol("$(param)_diff_app_minus_mapt")
        mean_col = Symbol("$(param)_mean")
        df[!, diff_col] = df[!, app_col] .- df[!, mapt_col]
        df[!, mean_col] = (df[!, app_col] .+ df[!, mapt_col]) ./ 2
    end
    sort!(df, :condition_pair_rank)
    return df
end

function pearsonr(x, y)
    length(x) < 3 && return NaN
    sx = std(x)
    sy = std(y)
    (sx == 0 || sy == 0) && return NaN
    return cor(x, y)
end

function paired_stats(df::DataFrame, protein::AbstractString, param::AbstractString, filter_name::AbstractString, mask)
    sub = df[mask, :]
    app = sub[!, Symbol("$(param)_app")]
    mapt = sub[!, Symbol("$(param)_mapt")]
    diff = app .- mapt
    n = length(diff)
    mean_diff = n == 0 ? NaN : mean(diff)
    sd_diff = n <= 1 ? NaN : std(diff)
    se_diff = n <= 1 ? NaN : sd_diff / sqrt(n)
    t_stat = (n <= 1 || se_diff == 0 || isnan(se_diff)) ? NaN : mean_diff / se_diff
    p_t = isnan(t_stat) ? NaN : 2 * ccdf(TDist(n - 1), abs(t_stat))
    return (
        protein = protein,
        parameter = param,
        filter = filter_name,
        n_regions = n,
        app_mean = n == 0 ? NaN : mean(app),
        mapt_mean = n == 0 ? NaN : mean(mapt),
        mean_diff_app_minus_mapt = mean_diff,
        median_diff_app_minus_mapt = n == 0 ? NaN : median(diff),
        sd_diff = sd_diff,
        paired_t = t_stat,
        paired_t_p = p_t,
        pearson_app_mapt = pearsonr(app, mapt),
        frac_app_greater = n == 0 ? NaN : mean(diff .> 0),
    )
end

function summary_stats(df::DataFrame, protein::AbstractString)
    rows = []
    for param in ["alpha", "beta", "gamma"]
        rhat_app = Symbol("$(param)_rhat_app")
        rhat_mapt = Symbol("$(param)_rhat_mapt")
        masks = Dict(
            "all" => trues(nrow(df)),
            "active_any" => df.active_any,
            "active_and_param_rhat_le_1_05" => df.active_any .& coalesce.(df[!, rhat_app] .<= 1.05, false) .& coalesce.(df[!, rhat_mapt] .<= 1.05, false),
        )
        for (filter_name, mask) in masks
            push!(rows, paired_stats(df, protein, param, filter_name, mask))
        end
    end
    return DataFrame(rows)
end

function scatter_identity!(plt, x, y)
    finite = isfinite.(x) .& isfinite.(y)
    if any(finite)
        minxy = min(minimum(x[finite]), minimum(y[finite]))
        maxxy = max(maximum(x[finite]), maximum(y[finite]))
        plot!(plt, [minxy, maxxy], [minxy, maxxy]; color = :black, linestyle = :dash, linewidth = 1.8, label = false)
    end
end

function parameter_scatter_panels(df::DataFrame, protein::AbstractString, out_path::AbstractString)
    params = ["alpha", "beta", "gamma"]
    subplots = Plots.Plot[]
    active = df.active_any
    for param in params
        app = df[!, Symbol("$(param)_app")]
        mapt = df[!, Symbol("$(param)_mapt")]
        r = pearsonr(app[active], mapt[active])
        plt = scatter(
            mapt[.!active],
            app[.!active];
            xlabel = "MAPT",
            ylabel = "APP",
            title = "$(param) (active r=$(round(r; digits = 2)))",
            label = "inactive",
            color = :gray70,
            markersize = 3,
            alpha = 0.55,
            markerstrokewidth = 0,
        )
        scatter!(
            plt,
            mapt[active],
            app[active];
            label = "active",
            color = RGB(0 / 255, 71 / 255, 171 / 255),
            markersize = 4,
            alpha = 0.75,
            markerstrokecolor = :white,
            markerstrokewidth = 0.4,
        )
        scatter_identity!(plt, mapt, app)
        push!(subplots, plt)
    end
    panel = plot(subplots...; layout = (1, 3), size = (1350, 420), plot_title = "$(uppercase(protein)) REGION-RF: APP vs MAPT")
    savefig(panel, out_path * ".pdf")
    savefig(panel, out_path * ".png")
end

function difference_by_rank_panels(df::DataFrame, protein::AbstractString, out_path::AbstractString)
    params = ["alpha", "beta", "gamma"]
    subplots = Plots.Plot[]
    active = df.active_any
    for param in params
        diff = df[!, Symbol("$(param)_diff_app_minus_mapt")]
        plt = scatter(
            df.condition_pair_rank[.!active],
            diff[.!active];
            xlabel = "Pathology rank across conditions",
            ylabel = "APP - MAPT",
            title = param,
            label = "inactive",
            color = :gray70,
            markersize = 3,
            alpha = 0.55,
            markerstrokewidth = 0,
        )
        scatter!(
            plt,
            df.condition_pair_rank[active],
            diff[active];
            label = "active",
            color = RGB(196 / 255, 54 / 255, 22 / 255),
            markersize = 4,
            alpha = 0.75,
            markerstrokecolor = :white,
            markerstrokewidth = 0.4,
        )
        hline!(plt, [0.0]; color = :black, linestyle = :dash, linewidth = 1.6, label = false)
        push!(subplots, plt)
    end
    panel = plot(subplots...; layout = (1, 3), size = (1350, 420), plot_title = "$(uppercase(protein)) REGION-RF: APP - MAPT by Pathology Rank")
    savefig(panel, out_path * ".pdf")
    savefig(panel, out_path * ".png")
end

function difference_hist_panels(df::DataFrame, protein::AbstractString, out_path::AbstractString)
    params = ["alpha", "beta", "gamma"]
    subplots = Plots.Plot[]
    for param in params
        diff = df[df.active_any, Symbol("$(param)_diff_app_minus_mapt")]
        plt = histogram(
            diff;
            xlabel = "APP - MAPT",
            ylabel = "Active regions",
            title = param,
            bins = 35,
            legend = false,
            color = RGB(0 / 255, 71 / 255, 171 / 255),
            alpha = 0.75,
        )
        vline!(plt, [0.0]; color = :black, linestyle = :dash, linewidth = 1.6)
        vline!(plt, [median(diff)]; color = RGB(196 / 255, 54 / 255, 22 / 255), linewidth = 2.0)
        push!(subplots, plt)
    end
    panel = plot(subplots...; layout = (1, 3), size = (1350, 420), plot_title = "$(uppercase(protein)) REGION-RF active-region parameter shifts")
    savefig(panel, out_path * ".pdf")
    savefig(panel, out_path * ".png")
end

function make_plots(df::DataFrame, protein::AbstractString, figure_dir::AbstractString)
    mkpath(figure_dir)
    parameter_scatter_panels(df, protein, joinpath(figure_dir, "$(protein)_app_vs_mapt_parameters"))
    difference_by_rank_panels(df, protein, joinpath(figure_dir, "$(protein)_app_minus_mapt_by_rank"))
    difference_hist_panels(df, protein, joinpath(figure_dir, "$(protein)_active_difference_histograms"))
end

function main()
    project_root = abspath(get_arg("--project-root", dirname(dirname(@__DIR__))))
    out_dir = get_arg("--out-dir", joinpath(project_root, "paper-copath", "results", "region_rf_condition_comparison"))
    figure_dir = get_arg("--figure-dir", joinpath(project_root, "paper-copath", "figures", "region_rf_condition_comparison"))
    mkpath(out_dir)
    mkpath(figure_dir)

    summaries = DataFrame()
    for protein in ["syn", "tau"]
        df = comparison_table(project_root, protein)
        CSV.write(joinpath(out_dir, "$(protein)_app_vs_mapt_region_parameters.csv"), df)
        stats = summary_stats(df, protein)
        append!(summaries, stats; cols = :union)
        make_plots(df, protein, figure_dir)
    end
    CSV.write(joinpath(out_dir, "app_vs_mapt_parameter_shift_summary.csv"), summaries)
    println(out_dir)
end

main()
