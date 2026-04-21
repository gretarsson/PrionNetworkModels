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
        return (samples = samples, parameter_names = parameter_names)
    end
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
    else
        error("Run bundle is missing predictions_train.csv: $run_dir")
    end

    return outdir
end
