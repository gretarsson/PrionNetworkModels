function load_run_matrix(path::AbstractString)
    df = CSV.read(path, DataFrame)
    labels = String.(df[:, 1])
    timepoint_names = String.(names(df)[2:end])
    timepoints = parse.(Float64, timepoint_names)
    matrix = Matrix{Float64}(df[:, 2:end])
    return (labels = labels, timepoints = timepoints, values = matrix)
end

function load_posterior_draws(path::AbstractString)
    h5open(path, "r") do h5
        samples = Matrix{Float64}(read(h5["chains/samples"]))
        parameter_names = String.(read(h5["chains/parameter_names"]))
        chain_ids = haskey(h5, "chains/chain_ids") ? Int.(read(h5["chains/chain_ids"])) : ones(Int, size(samples, 1))
        return (samples = samples, parameter_names = parameter_names, chain_ids = chain_ids)
    end
end

function posterior_chains(posterior)
    chain_ids = Int.(posterior.chain_ids)
    unique_ids = sort(unique(chain_ids))
    length(unique_ids) >= 2 || error("Need at least two chains for convergence diagnostics")

    sample_counts = [count(==(chain_id), chain_ids) for chain_id in unique_ids]
    n_samples = minimum(sample_counts)
    n_params = size(posterior.samples, 2)
    arr = Array{Float64}(undef, n_samples, n_params, length(unique_ids))

    for (j, chain_id) in enumerate(unique_ids)
        idxs = findall(==(chain_id), chain_ids)
        arr[:, :, j] = posterior.samples[idxs[1:n_samples], :]
    end

    return Chains(arr, Symbol.(posterior.parameter_names))
end

function chains_long_df(chain::Chains)
    arr = Array(chain)
    n_iter, n_param, n_chain = size(arr)
    params = String.(names(chain, :parameters))

    iterations = Int[]
    chain_ids = Int[]
    parameter_names = String[]
    values = Float64[]
    sizehint!(iterations, n_iter * n_param * n_chain)
    sizehint!(chain_ids, n_iter * n_param * n_chain)
    sizehint!(parameter_names, n_iter * n_param * n_chain)
    sizehint!(values, n_iter * n_param * n_chain)

    for c in 1:n_chain
        for p in 1:n_param
            pname = params[p]
            for i in 1:n_iter
                push!(iterations, i)
                push!(chain_ids, c)
                push!(parameter_names, pname)
                push!(values, arr[i, p, c])
            end
        end
    end

    return DataFrame(iteration = iterations, chain = chain_ids, parameter = parameter_names, value = values)
end

is_local_param(name::String) = startswith(name, "beta[") || startswith(name, "gamma[")

function split_local_params(rhats::Dict{String,Float64})
    beta = Dict{String,Float64}()
    gamma = Dict{String,Float64}()
    for (name, val) in rhats
        if startswith(name, "beta[")
            beta[name] = val
        elseif startswith(name, "gamma[")
            gamma[name] = val
        end
    end
    return beta, gamma
end

function compute_rhat_semantic(chain::Chains)
    rhat_obj = MCMCChains.MCMCDiagnosticTools.rhat(chain)
    raw_names = String.(rhat_obj.nt.parameters)
    raw_vals = rhat_obj.nt.rhat

    semantic_rhat = Dict{String,Float64}()
    for (name, val) in zip(raw_names, raw_vals)
        if name == "lp" || name == "lp__"
            continue
        end
        semantic_rhat[name] = val
    end
    return semantic_rhat
end

function compute_summary_stats(chain::Chains)
    return DataFrame(MCMCChains.summarystats(chain))
end

function top_problem_parameters(summary::DataFrame; top_k_local::Int=12)
    summary = copy(summary)
    summary.parameter = String.(summary.parameters)
    summary.family = map(summary.parameter) do name
        if startswith(name, "beta[")
            "beta"
        elseif startswith(name, "gamma[")
            "gamma"
        else
            "global"
        end
    end

    globals = sort(summary[summary.family .== "global", :], :rhat, rev = true)
    locals = sort(summary[summary.family .!= "global", :], :rhat, rev = true)
    selected = vcat(globals, first(locals, min(top_k_local, nrow(locals))))
    return unique(selected[:, [:parameter, :family, :rhat, :ess_bulk, :ess_tail]])
end

function plot_trace(output_path::AbstractString, chain::Chains, parameter::AbstractString)
    symbol_name = Symbol(parameter)
    arr = Array(chain[symbol_name])
    ndims(arr) == 2 || error("Unexpected chain shape for $parameter")
    n_iter, n_chain = size(arr)

    plt = plot(
        xlabel = "Iteration",
        ylabel = parameter,
        title = "Trace: $parameter",
        legend = :topright,
        size = (900, 500),
    )

    colors = [
        RGB(0 / 255, 71 / 255, 171 / 255),
        RGB(196 / 255, 54 / 255, 22 / 255),
        RGB(0 / 255, 136 / 255, 55 / 255),
        RGB(123 / 255, 31 / 255, 162 / 255),
    ]
    for c in 1:n_chain
        plot!(
            plt,
            1:n_iter,
            arr[:, c];
            label = "Chain $c",
            linewidth = 1.7,
            alpha = 0.9,
            color = colors[mod1(c, length(colors))],
        )
    end

    savefig(plt, output_path)
    return output_path
end

function chain_parameter_means(chain::Chains, parameters::Vector{String})
    rows = DataFrame(parameter = String[], chain = Int[], mean = Float64[], sd = Float64[])
    for parameter in parameters
        arr = Array(chain[Symbol(parameter)])
        n_iter, n_chain = size(arr)
        for c in 1:n_chain
            push!(rows, (parameter, c, mean(arr[:, c]), std(arr[:, c])))
        end
    end
    return rows
end

function chain_fit_metrics(run_dir::AbstractString)
    spec = resolve_bundle_spec_paths(load_run_spec(joinpath(run_dir, "spec.toml")), run_dir)
    transport = build_transport_operator(spec.data.network; transport = spec.model.transport)
    pathology = process_pathology(spec.data.observations; network_csv = spec.data.network)
    posterior = load_posterior_draws(joinpath(run_dir, "posterior.h5"))

    param_names = posterior.parameter_names
    samples = posterior.samples
    chain_ids = posterior.chain_ids
    name_to_idx = Dict(name => idx for (idx, name) in enumerate(param_names))
    n_regions = length(transport.labels)
    parameter_order = trajectory_parameter_names(spec.model.name, n_regions)
    parameter_indices = [get(name_to_idx, name, 0) for name in parameter_order]
    any(==(0), parameter_indices) && error("Posterior draws are missing required parameters for $(spec.model.name)")

    seed_index = get(name_to_idx, "seed", 0)
    sigma_index = get(name_to_idx, "sigma", 0)
    (seed_index == 0 || sigma_index == 0) && error("Posterior draws are missing seed or sigma parameters")

    ignore_regions = spec.inference.ignore_seed ? spec.seeding.seed_indices : Int[]
    obs_all = finite_observations(pathology.data; mean_data = spec.inference.mean_data, ignore_regions = ignore_regions)
    seed_data = pathology.data[spec.seeding.seed_indices, :, :]
    obs_seed = finite_observations(seed_data; mean_data = spec.inference.mean_data, ignore_regions = Int[])

    rows = DataFrame(
        chain = Int[],
        loglik_all = Float64[],
        loglik_seed = Float64[],
        rmse_all = Float64[],
        mae_all = Float64[],
        sigma = Float64[],
        seed = Float64[],
    )

    for chain_id in sort(unique(chain_ids))
        idxs = findall(==(chain_id), chain_ids)
        chain_mean = vec(mean(samples[idxs, :]; dims = 1))
        params = collect(chain_mean[parameter_indices])
        seed_value = chain_mean[seed_index]
        sigma = chain_mean[sigma_index]

        predicted = simulate_trajectory(spec, transport.L, transport.labels, pathology.timepoints, params; seed_value = seed_value)
        predicted_obs = predicted[1:n_regions, :]

        if spec.inference.mean_data
            observed_matrix_all = observed_matrix(pathology.data)
            pred_all = vec(predicted_obs[obs_all.mask])
            obs_vec_all = obs_all.values

            observed_matrix_seed = observed_matrix(seed_data)
            pred_seed = vec(predicted_obs[spec.seeding.seed_indices, :][obs_seed.mask])
            obs_vec_seed = obs_seed.values
        else
            pred_rep_all = cat([predicted_obs for _ in 1:obs_all.n_samples]..., dims = 3)
            pred_all = vec(pred_rep_all)[obs_all.nonmissing]
            obs_vec_all = obs_all.values

            pred_seed_regions = predicted_obs[spec.seeding.seed_indices, :]
            pred_rep_seed = cat([pred_seed_regions for _ in 1:obs_seed.n_samples]..., dims = 3)
            pred_seed = vec(pred_rep_seed)[obs_seed.nonmissing]
            obs_vec_seed = obs_seed.values
        end

        residuals = obs_vec_all .- pred_all
        loglik_all = sum(logpdf.(Normal.(pred_all, sigma), obs_vec_all))
        loglik_seed = sum(logpdf.(Normal.(pred_seed, sigma), obs_vec_seed))
        rmse_all = sqrt(mean(residuals .^ 2))
        mae_all = mean(abs.(residuals))

        push!(rows, (chain_id, loglik_all, loglik_seed, rmse_all, mae_all, sigma, seed_value))
    end

    return rows
end

function plot_chain_loglik(output_path::AbstractString, metrics::DataFrame)
    nrow(metrics) == 0 && return nothing

    chain_labels = string.(metrics.chain)
    xs = collect(1:nrow(metrics))

    plt_all = scatter(
        xs,
        metrics.loglik_all;
        xlabel = "Chain",
        ylabel = "Log-likelihood",
        title = "All Observations",
        markersize = 8,
        alpha = 0.9,
        color = RGB(0 / 255, 71 / 255, 171 / 255),
        legend = false,
        xticks = (xs, chain_labels),
    )

    plt_seed = scatter(
        xs,
        metrics.loglik_seed;
        xlabel = "Chain",
        ylabel = "Log-likelihood",
        title = "Seed Region Observations",
        markersize = 8,
        alpha = 0.9,
        color = RGB(196 / 255, 54 / 255, 22 / 255),
        legend = false,
        xticks = (xs, chain_labels),
    )

    plt = plot(plt_all, plt_seed; layout = (2, 1), size = (900, 800))
    savefig(plt, output_path)
    return output_path
end

function plot_seed_region_chain_comparison(output_path::AbstractString, run_dir::AbstractString)
    spec = resolve_bundle_spec_paths(load_run_spec(joinpath(run_dir, "spec.toml")), run_dir)
    transport = build_transport_operator(spec.data.network; transport = spec.model.transport)
    pathology = process_pathology(spec.data.observations; network_csv = spec.data.network)
    posterior = load_posterior_draws(joinpath(run_dir, "posterior.h5"))

    param_names = posterior.parameter_names
    samples = posterior.samples
    chain_ids = posterior.chain_ids
    name_to_idx = Dict(name => idx for (idx, name) in enumerate(param_names))
    n_regions = length(transport.labels)
    parameter_order = trajectory_parameter_names(spec.model.name, n_regions)
    parameter_indices = [get(name_to_idx, name, 0) for name in parameter_order]
    any(==(0), parameter_indices) && error("Posterior draws are missing required parameters for $(spec.model.name)")

    seed_index = get(name_to_idx, "seed", 0)
    seed_index == 0 && error("Posterior draws are missing the seed parameter")

    summary = summarize_over_replicates(pathology.data)
    seed_regions = spec.seeding.seed_indices
    seed_labels = pathology.labels[seed_regions]
    timepoints = pathology.timepoints

    colors = [
        RGB(0 / 255, 71 / 255, 171 / 255),
        RGB(196 / 255, 54 / 255, 22 / 255),
        RGB(0 / 255, 136 / 255, 55 / 255),
        RGB(123 / 255, 31 / 255, 162 / 255),
    ]

    subplots = Plots.Plot[]

    for (panel_idx, region_idx) in enumerate(seed_regions)
        subplot = plot(
            xlabel = "Time",
            ylabel = "Pathology",
            title = "Seed Region: $(seed_labels[panel_idx])",
            legend = :topright,
            ylims = (0.0, max(maximum(skipmissing(summary.mean[region_idx, :])), 1e-3) * 1.15),
        )

        scatter!(
            subplot,
            timepoints,
            summary.mean[region_idx, :];
            yerror = summary.se[region_idx, :],
            markersize = 5,
            color = :black,
            label = "Observed mean ± SE",
        )

        for chain_id in sort(unique(chain_ids))
            idxs = findall(==(chain_id), chain_ids)
            chain_mean = vec(mean(samples[idxs, :]; dims = 1))
            params = collect(chain_mean[parameter_indices])
            seed_value = chain_mean[seed_index]
            predicted = simulate_trajectory(spec, transport.L, transport.labels, timepoints, params; seed_value = seed_value)
            pred_obs = predicted[1:n_regions, :]
            plot!(
                subplot,
                timepoints,
                pred_obs[region_idx, :];
                label = "Chain $chain_id",
                linewidth = 2.0,
                alpha = 0.9,
                color = colors[mod1(chain_id, length(colors))],
            )
        end

        push!(subplots, subplot)
    end

    plt = plot(subplots...; layout = (length(seed_regions), 1), size = (900, 320 * length(seed_regions)))
    savefig(plt, output_path)
    return output_path
end

function plot_rhat_scatter(output_path::AbstractString, rhats::Dict{String,Float64}; title::AbstractString="")
    isempty(rhats) && return nothing

    selected = collect(keys(rhats))
    sort!(selected)
    ys = [rhats[name] for name in selected]
    xs = 1:length(selected)

    plt = scatter(
        xs,
        ys;
        xlabel = "Parameter index",
        ylabel = "Rhat",
        title = String(title),
        markersize = 6,
        alpha = 0.85,
        color = RGB(0 / 255, 71 / 255, 171 / 255),
        legend = false,
        size = (900, 500),
    )

    for (val, col, lw) in ((1.00, :gray60, 1.0), (1.01, :green4, 1.5), (1.05, :orange, 2.0), (1.10, :red, 2.0))
        hline!(plt, [val]; color = col, lw = lw, ls = :dash)
    end

    if length(selected) <= 25
        xticks!(plt, xs, selected)
        plot!(plt; xrotation = 45)
    end

    savefig(plt, output_path)
    return output_path
end

function diagnostics_plots(run_dir::AbstractString, output_dir::AbstractString)
    mkpath(output_dir)

    posterior = load_posterior_draws(joinpath(run_dir, "posterior.h5"))
    chain = posterior_chains(posterior)
    summary = compute_summary_stats(chain)
    rhats = compute_rhat_semantic(chain)

    global_rhats = Dict(name => val for (name, val) in rhats if !is_local_param(name))
    beta_rhats, gamma_rhats = split_local_params(rhats)

    plot_rhat_scatter(joinpath(output_dir, "global_rhat.pdf"), global_rhats; title = "Global Parameters")
    isempty(beta_rhats) || plot_rhat_scatter(joinpath(output_dir, "beta_rhat.pdf"), beta_rhats; title = "Beta Parameters")
    isempty(gamma_rhats) || plot_rhat_scatter(joinpath(output_dir, "gamma_rhat.pdf"), gamma_rhats; title = "Gamma Parameters")

    names = String[]
    families = String[]
    values = Float64[]
    for (name, val) in sort(collect(global_rhats); by = first)
        push!(names, name); push!(families, "global"); push!(values, val)
    end
    for (name, val) in sort(collect(beta_rhats); by = first)
        push!(names, name); push!(families, "beta"); push!(values, val)
    end
    for (name, val) in sort(collect(gamma_rhats); by = first)
        push!(names, name); push!(families, "gamma"); push!(values, val)
    end
    CSV.write(joinpath(output_dir, "rhat_summary.csv"), DataFrame(parameter = names, family = families, rhat = values))

    CSV.write(joinpath(output_dir, "summary_stats.csv"), summary)

    flagged = top_problem_parameters(summary)
    CSV.write(joinpath(output_dir, "flagged_parameters.csv"), flagged)

    trace_dir = joinpath(output_dir, "trace")
    mkpath(trace_dir)
    for parameter in flagged.parameter
        safe_name = replace(parameter, "[" => "_", "]" => "", "/" => "_")
        plot_trace(joinpath(trace_dir, "trace_$(safe_name).pdf"), chain, parameter)
    end

    chain_means = chain_parameter_means(chain, collect(flagged.parameter))
    CSV.write(joinpath(output_dir, "flagged_chain_means.csv"), chain_means)

    fit_metrics = chain_fit_metrics(run_dir)
    CSV.write(joinpath(output_dir, "chain_fit_metrics.csv"), fit_metrics)
    plot_chain_loglik(joinpath(output_dir, "chain_loglik_comparison.pdf"), fit_metrics)
    plot_seed_region_chain_comparison(joinpath(output_dir, "seed_region_chain_comparison.pdf"), run_dir)

    return output_dir
end

function load_observation_summary(path_csv::AbstractString; network_csv::Union{Nothing,String}=nothing)
    pathology = process_pathology(path_csv; network_csv = network_csv)
    summary = summarize_over_replicates(pathology.data)
    return (
        labels = pathology.labels,
        timepoints = pathology.timepoints,
        mean = summary.mean,
        sd = summary.sd,
        se = summary.se,
        n = summary.n,
    )
end

function trajectory_parameter_names(model_name::AbstractString, n_regions::Integer)
    if model_name == "DIFF"
        return ["rho"]
    elseif model_name == "DIFF-R"
        return vcat(["rho", "alpha"], ["beta[$i]" for i in 1:n_regions])
    elseif model_name == "DIFF-RF"
        return vcat(["rho", "alpha"], ["beta[$i]" for i in 1:n_regions], ["gamma[$i]" for i in 1:n_regions])
    else
        error("Unknown model name: $model_name")
    end
end

function posterior_mean_retrodiction(run_dir::AbstractString, obs; n_dense::Int=300)
    spec = resolve_bundle_spec_paths(load_run_spec(joinpath(run_dir, "spec.toml")), run_dir)
    transport = build_transport_operator(spec.data.network; transport = spec.model.transport)
    posterior = load_posterior_draws(joinpath(run_dir, "posterior.h5"))
    size(posterior.samples, 1) == 0 && error("No posterior draws found in $(joinpath(run_dir, "posterior.h5"))")

    name_to_idx = Dict(name => idx for (idx, name) in enumerate(posterior.parameter_names))
    parameter_order = trajectory_parameter_names(spec.model.name, length(obs.labels))
    parameter_indices = [get(name_to_idx, name, 0) for name in parameter_order]
    any(==(0), parameter_indices) && error("Posterior draws are missing required parameters for $(spec.model.name)")

    seed_index = get(name_to_idx, "seed", 0)
    seed_index == 0 && error("Posterior draws are missing the seed parameter")

    sigma_index = get(name_to_idx, "sigma", 0)
    sigma_index == 0 && error("Posterior draws are missing the sigma parameter")

    mean_draw = vec(mean(posterior.samples; dims = 1))
    params = collect(mean_draw[parameter_indices])
    seed_value = mean_draw[seed_index]
    sigma = mean_draw[sigma_index]

    max_time = maximum(obs.timepoints)
    dense_timepoints = collect(range(0.0, max_time; length = n_dense))
    predicted = simulate_trajectory(spec, transport.L, transport.labels, dense_timepoints, params; seed_value = seed_value)
    mean_path = predicted[1:length(obs.labels), :]

    inner_z = quantile(Normal(), 0.75)
    outer_z = quantile(Normal(), 0.95)

    lower50 = max.(mean_path .- inner_z .* sigma, 0.0)
    upper50 = mean_path .+ inner_z .* sigma
    lower90 = max.(mean_path .- outer_z .* sigma, 0.0)
    upper90 = mean_path .+ outer_z .* sigma

    return (
        timepoints = dense_timepoints,
        mean = mean_path,
        lower50 = lower50,
        upper50 = upper50,
        lower90 = lower90,
        upper90 = upper90,
        sigma = sigma,
    )
end

function predicted_observed_plot(obs, predicted_path::AbstractString, output_path::AbstractString)
    pred = load_run_matrix(predicted_path)
    x = vec(obs.mean)
    y = vec(pred.values)
    finite = isfinite.(x) .& isfinite.(y)
    x = x[finite]
    y = y[finite]
    minxy = min(minimum(x), minimum(y))
    maxxy = max(maximum(x), maximum(y))

    plt = scatter(
        x,
        y;
        xlabel = "Observed",
        ylabel = "Predicted",
        title = "Predicted vs Observed",
        legend = false,
        alpha = 0.7,
        markersize = 4,
        color = RGB(0 / 255, 71 / 255, 171 / 255),
        markerstrokecolor = :white,
        markerstrokewidth = 0.5,
    )
    plot!(plt, [minxy, maxxy], [minxy, maxxy]; color = :black, linestyle = :dash, linewidth = 2)
    savefig(plt, output_path)
    return output_path
end

function retrodiction_plots(obs, output_dir::AbstractString; run_dir::Union{Nothing,String}=nothing)
    mkpath(output_dir)

    smooth = nothing
    if !isnothing(run_dir) && isfile(joinpath(run_dir, "posterior.h5")) && isfile(joinpath(run_dir, "spec.toml"))
        smooth = posterior_mean_retrodiction(String(run_dir), obs)
    end

    global_ymax = maximum(skipmissing(vec(obs.mean)))
    if !isnothing(smooth)
        global_ymax = max(global_ymax, maximum(smooth.upper90))
    end
    global_ymax = max(global_ymax, 1e-3) * 1.05

    observed_color = RGB(0 / 255, 71 / 255, 171 / 255)
    for i in eachindex(obs.labels)
        plt = plot(
            xlabel = "Time",
            ylabel = "Pathology",
            title = obs.labels[i],
            legend = :topleft,
            ylims = (0.0, global_ymax),
            linewidth = 3,
            size = (720, 420),
        )

        if !isnothing(smooth)
            plot!(
                plt,
                smooth.timepoints,
                smooth.mean[i, :];
                ribbon = (smooth.mean[i, :] .- smooth.lower90[i, :], smooth.upper90[i, :] .- smooth.mean[i, :]),
                fillalpha = 0.18,
                fillcolor = :gray70,
                linealpha = 0.0,
                label = "90% noise band",
            )
            plot!(
                plt,
                smooth.timepoints,
                smooth.mean[i, :];
                ribbon = (smooth.mean[i, :] .- smooth.lower50[i, :], smooth.upper50[i, :] .- smooth.mean[i, :]),
                fillalpha = 0.28,
                fillcolor = :gray45,
                linealpha = 0.0,
                label = "50% noise band",
            )
            plot!(
                plt,
                smooth.timepoints,
                smooth.mean[i, :];
                color = :black,
                linewidth = 3,
                label = "Posterior mean",
            )
        end

        scatter!(
            plt,
            obs.timepoints,
            obs.mean[i, :];
            label = "Observed mean",
            color = observed_color,
            markersize = 5,
            markerstrokecolor = :white,
            markerstrokewidth = 0.7,
        )

        savefig(plt, joinpath(output_dir, "retrodiction_$(obs.labels[i]).pdf"))
    end
    return output_dir
end

function plot_run_bundle(run_dir::AbstractString; output_dir::Union{Nothing,String}=nothing)
    outdir = isnothing(output_dir) ? joinpath(run_dir, "plots") : String(output_dir)
    mkpath(outdir)

    spec = resolve_bundle_spec_paths(load_run_spec(joinpath(run_dir, "spec.toml")), run_dir)
    predictions_train = joinpath(run_dir, "predictions_train.csv")
    obs = load_observation_summary(spec.data.observations; network_csv = spec.data.network)

    if isfile(predictions_train)
        predicted_observed_plot(
            obs,
            predictions_train,
            joinpath(outdir, "predicted_vs_observed.pdf"),
        )
        retrodiction_plots(
            obs,
            joinpath(outdir, "retrodiction");
            run_dir = run_dir,
        )
        if isfile(joinpath(run_dir, "posterior.h5"))
            diagnostics_plots(run_dir, joinpath(outdir, "diagnostics"))
        end
    else
        error("Run bundle is missing predictions_train.csv: $run_dir")
    end

    return outdir
end
