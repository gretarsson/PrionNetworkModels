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

function rank_average(values::AbstractVector)
    order = sortperm(values)
    ranks = zeros(Float64, length(values))
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && values[order[j + 1]] == values[order[i]]
            j += 1
        end
        avg_rank = (i + j) / 2
        for k in i:j
            ranks[order[k]] = avg_rank
        end
        i = j + 1
    end
    return ranks
end

function spearmanr(x, y)
    length(x) < 3 && return NaN
    return pearsonr(rank_average(x), rank_average(y))
end

function correlation_pvalue(r::Real, n::Integer)
    if n < 3 || !isfinite(r) || abs(r) >= 1
        return NaN
    end
    t = r * sqrt((n - 2) / (1 - r^2))
    return 2 * ccdf(TDist(n - 2), abs(t))
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

function maybe_float(value)
    ismissing(value) && return missing
    value isa Number && return Float64(value)
    parsed = tryparse(Float64, String(value))
    return isnothing(parsed) ? missing : parsed
end

function region_means_from_wide(path::AbstractString)
    df = CSV.read(path, DataFrame)
    means = Dict{String,Float64}()
    for col in names(df)[3:end]
        values = collect(skipmissing(maybe_float.(df[!, col])))
        values = values[isfinite.(values)]
        isempty(values) && continue
        means[String(col)] = mean(values)
    end
    return means
end

function add_amyloid_columns!(df::DataFrame, project_root::AbstractString, protein::AbstractString)
    data_dir = joinpath(project_root, "paper-copath", "data")
    # Available amyloid tables are from APP/MAPTApp KI mice only, with injected
    # treatment and non-injected controls. There is no separate MAPT amyloid
    # condition in the processed files, so we compare treatment minus control.
    treatment = protein == "syn" ? "mpff" : "adphf"
    ab40_treatment = region_means_from_wide(joinpath(data_dir, "ab40_pathology_$(treatment).csv"))
    ab42_treatment = region_means_from_wide(joinpath(data_dir, "ab42_pathology_$(treatment).csv"))
    ab40_control = region_means_from_wide(joinpath(data_dir, "ab40_pathology_control.csv"))
    ab42_control = region_means_from_wide(joinpath(data_dir, "ab42_pathology_control.csv"))

    df.abeta_treatment = fill(treatment, nrow(df))
    df.ab40_treatment_mean = [get(ab40_treatment, r, NaN) for r in df.region]
    df.ab42_treatment_mean = [get(ab42_treatment, r, NaN) for r in df.region]
    df.ab40_control_mean = [get(ab40_control, r, NaN) for r in df.region]
    df.ab42_control_mean = [get(ab42_control, r, NaN) for r in df.region]
    df.ab40_diff_treatment_minus_control = df.ab40_treatment_mean .- df.ab40_control_mean
    df.ab42_diff_treatment_minus_control = df.ab42_treatment_mean .- df.ab42_control_mean
    df.ab40_diff_log10p1 = log10.(1 .+ max.(df.ab40_treatment_mean, 0.0)) .- log10.(1 .+ max.(df.ab40_control_mean, 0.0))
    df.ab42_diff_log10p1 = log10.(1 .+ max.(df.ab42_treatment_mean, 0.0)) .- log10.(1 .+ max.(df.ab42_control_mean, 0.0))
    return df
end

function amyloid_correlation_stats(df::DataFrame, protein::AbstractString)
    rows = []
    for amyloid in ["ab40", "ab42"]
        xcol = Symbol("$(amyloid)_diff_log10p1")
        for param in ["alpha", "beta", "gamma"]
            ycol = Symbol("$(param)_diff_app_minus_mapt")
            rhat_app = Symbol("$(param)_rhat_app")
            rhat_mapt = Symbol("$(param)_rhat_mapt")
            filters = Dict(
                "all" => trues(nrow(df)),
                "active_any" => df.active_any,
                "active_and_param_rhat_le_1_05" => df.active_any .& coalesce.(df[!, rhat_app] .<= 1.05, false) .& coalesce.(df[!, rhat_mapt] .<= 1.05, false),
            )
            for (filter_name, filter_mask) in filters
                finite = filter_mask .& isfinite.(df[!, xcol]) .& isfinite.(df[!, ycol])
                x = df[finite, xcol]
                y = df[finite, ycol]
                pearson = pearsonr(x, y)
                spearman = spearmanr(x, y)
                push!(rows, (
                    protein = protein,
                    amyloid = amyloid,
                    parameter = param,
                    filter = filter_name,
                    n_regions = length(x),
                    pearson_r = pearson,
                    pearson_p = correlation_pvalue(pearson, length(x)),
                    spearman_r = spearman,
                    spearman_p = correlation_pvalue(spearman, length(x)),
                    amyloid_mean = isempty(x) ? NaN : mean(x),
                    parameter_shift_mean = isempty(y) ? NaN : mean(y),
                ))
            end
        end
    end
    return DataFrame(rows)
end

function scatter_with_fit!(plt, x, y)
    finite = isfinite.(x) .& isfinite.(y)
    count(finite) < 3 && return plt
    xf = x[finite]
    yf = y[finite]
    X = hcat(ones(length(xf)), xf)
    coef = X \ yf
    xs = collect(range(minimum(xf), maximum(xf); length = 200))
    plot!(plt, xs, coef[1] .+ coef[2] .* xs; color = :black, linewidth = 2.0, label = false)
    return plt
end

function amyloid_scatter_panels(df::DataFrame, protein::AbstractString, amyloid::AbstractString, figure_dir::AbstractString)
    xcol = Symbol("$(amyloid)_diff_log10p1")
    params = ["alpha", "beta", "gamma"]
    active = df.active_any
    subplots = Plots.Plot[]
    for param in params
        ycol = Symbol("$(param)_diff_app_minus_mapt")
        x = df[!, xcol]
        y = df[!, ycol]
        mask = active .& isfinite.(x) .& isfinite.(y)
        r = pearsonr(x[mask], y[mask])
        p = correlation_pvalue(r, count(mask))
        plt = scatter(
            x[.!active],
            y[.!active];
            xlabel = "$(uppercase(amyloid)) treatment - control log10(1 + burden)",
            ylabel = "$(param) APP - MAPT",
            title = "$(param): r=$(round(r; digits = 2)), p=$(round(p; sigdigits = 2))",
            label = "inactive",
            color = :gray70,
            markersize = 3,
            alpha = 0.45,
            markerstrokewidth = 0,
        )
        scatter!(
            plt,
            x[active],
            y[active];
            label = "active",
            color = RGB(0 / 255, 71 / 255, 171 / 255),
            markersize = 4,
            alpha = 0.75,
            markerstrokecolor = :white,
            markerstrokewidth = 0.4,
        )
        hline!(plt, [0.0]; color = :black, linestyle = :dash, linewidth = 1.4, label = false)
        scatter_with_fit!(plt, x[mask], y[mask])
        push!(subplots, plt)
    end
    panel = plot(
        subplots...;
        layout = (1, 3),
        size = (1400, 430),
        plot_title = "$(uppercase(protein)) REGION-RF APP - MAPT shifts vs $(uppercase(amyloid))",
    )
    out = joinpath(figure_dir, "$(protein)_$(amyloid)_amyloid_vs_parameter_shifts")
    savefig(panel, out * ".pdf")
    savefig(panel, out * ".png")
end

function amyloid_plots(df::DataFrame, protein::AbstractString, figure_dir::AbstractString)
    for amyloid in ["ab40", "ab42"]
        amyloid_scatter_panels(df, protein, amyloid, figure_dir)
    end
end

function main()
    project_root = abspath(get_arg("--project-root", dirname(dirname(@__DIR__))))
    out_dir = get_arg("--out-dir", joinpath(project_root, "paper-copath", "results", "region_rf_condition_comparison"))
    figure_dir = get_arg("--figure-dir", joinpath(project_root, "paper-copath", "figures", "region_rf_condition_comparison"))
    mkpath(out_dir)
    mkpath(figure_dir)

    summaries = DataFrame()
    amyloid_summaries = DataFrame()
    for protein in ["syn", "tau"]
        df = comparison_table(project_root, protein)
        add_amyloid_columns!(df, project_root, protein)
        CSV.write(joinpath(out_dir, "$(protein)_app_vs_mapt_region_parameters.csv"), df)
        stats = summary_stats(df, protein)
        append!(summaries, stats; cols = :union)
        append!(amyloid_summaries, amyloid_correlation_stats(df, protein); cols = :union)
        make_plots(df, protein, figure_dir)
        amyloid_plots(df, protein, figure_dir)
    end
    CSV.write(joinpath(out_dir, "app_vs_mapt_parameter_shift_summary.csv"), summaries)
    CSV.write(joinpath(out_dir, "amyloid_vs_parameter_shift_correlations.csv"), amyloid_summaries)
    println(out_dir)
end

main()
