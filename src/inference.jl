function default_priors(model_name::AbstractString, N::Integer)
    name = String(model_name)
    if name == "DIFF"
        return (
            rho = truncated(Normal(0.2, 0.1), 0.0, Inf),
            sigma = LogNormal(0, 1),
            seed = truncated(Normal(0.0, 0.1), 0.0, Inf),
        )
    elseif name == "DIFF-R"
        return (
            rho = truncated(Normal(0.0, 0.1), 0.0, Inf),
            alpha = truncated(Normal(0.0, 0.1), 0.0, Inf),
            beta = Normal(0.0, 1.0),
            sigma = LogNormal(0, 1),
            seed = truncated(Normal(0.0, 0.1), 0.0, Inf),
        )
    elseif name == "DIFF-RF"
        return (
            rho = truncated(Normal(0.0, 0.1), 0.0, Inf),
            alpha = truncated(Normal(0.0, 0.1), 0.0, Inf),
            beta = Normal(0.0, 1.0),
            gamma = truncated(Normal(0.0, 0.1), 0.0, Inf),
            sigma = LogNormal(0, 1),
            seed = truncated(Normal(0.0, 0.1), 0.0, Inf),
        )
    else
        error("Unknown model name: $model_name")
    end
end

function finite_observations(data::AbstractArray; mean_data::Bool=false, ignore_regions::Vector{Int}=Int[])
    working = copy(data)
    if !isempty(ignore_regions)
        if ndims(working) == 2
            working[ignore_regions, :] .= missing
        else
            working[ignore_regions, :, :] .= missing
        end
    end

    if mean_data || ndims(working) == 2
        obs = observed_matrix(working)
        mask = isfinite.(obs)
        vals = Float64.(obs[mask])
        return (values = vals, mask = mask, n_samples = 1, mean_data = true)
    end

    vec_data = vec(working)
    nonmissing = findall(vec_data .!== missing)
    vals = Float64.(vec_data[nonmissing])
    return (
        values = vals,
        nonmissing = nonmissing,
        n_samples = size(working, 3),
        mean_data = false,
    )
end

function parameter_names_for_model(model_name::AbstractString, N::Integer)
    name = String(model_name)
    if name == "DIFF"
        return ["rho", "sigma", "seed"]
    elseif name == "DIFF-R"
        return vcat(["rho", "alpha"], ["beta[$i]" for i in 1:N], ["sigma", "seed"])
    elseif name == "DIFF-RF"
        return vcat(["rho", "alpha"], ["beta[$i]" for i in 1:N], ["gamma[$i]" for i in 1:N], ["sigma", "seed"])
    else
        error("Unknown model name: $model_name")
    end
end

@model function smoke_model(spec::RunSpec, prob::ODEProblem, N::Integer, timepoints::Vector{Float64}, obs_vals::Vector{Float64}, obs_info)
    priors = default_priors(spec.model.name, N)

    if spec.model.name == "DIFF"
        rho ~ priors.rho
        sigma ~ priors.sigma
        seed ~ priors.seed
        params = [rho]
    elseif spec.model.name == "DIFF-R"
        rho ~ priors.rho
        alpha ~ priors.alpha
        beta ~ filldist(priors.beta, N)
        sigma ~ priors.sigma
        seed ~ priors.seed
        params = vcat([rho, alpha], beta)
    elseif spec.model.name == "DIFF-RF"
        rho ~ priors.rho
        alpha ~ priors.alpha
        beta ~ filldist(priors.beta, N)
        gamma ~ filldist(priors.gamma, N)
        sigma ~ priors.sigma
        seed ~ priors.seed
        params = vcat([rho, alpha], beta, gamma)
    else
        error("Unsupported model in smoke_model: $(spec.model.name)")
    end

    predicted = solve(
        prob,
        Tsit5();
        u0 = initial_conditions_for_spec(spec, N; seed_value = seed),
        p = collect(params),
        saveat = timepoints,
        sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP(true)),
        abstol = 1e-6,
        reltol = 1e-3,
        maxiters = 6000,
    )
    pred = Array(predicted)[1:N, :]
    pred_vec = if obs_info.mean_data
        vec(pred[obs_info.mask])
    else
        pred_rep = cat([pred for _ in 1:obs_info.n_samples]..., dims = 3)
        vec(pred_rep)[obs_info.nonmissing]
    end
    obs_vals ~ MvNormal(pred_vec, sigma^2 * I)
end

function posterior_summary_table(chain::Chains, names_order::Vector{String})
    out_names = String[]
    means = Float64[]
    sds = Float64[]
    for name in names_order
        samples = vec(Array(chain[Symbol(name)]))
        push!(out_names, name)
        push!(means, mean(samples))
        push!(sds, std(samples))
    end
    return DataFrame(parameter = out_names, mean = means, sd = sds)
end

function posterior_mean_parameter_vector(chain::Chains, model_name::AbstractString, N::Integer)
    name = String(model_name)
    if name == "DIFF"
        return [mean(vec(Array(chain[:rho])))]
    elseif name == "DIFF-R"
        beta = [mean(vec(Array(chain[Symbol("beta[$i]")]))) for i in 1:N]
        return vcat([mean(vec(Array(chain[:rho]))), mean(vec(Array(chain[:alpha])))], beta)
    elseif name == "DIFF-RF"
        beta = [mean(vec(Array(chain[Symbol("beta[$i]")]))) for i in 1:N]
        gamma = [mean(vec(Array(chain[Symbol("gamma[$i]")]))) for i in 1:N]
        return vcat([mean(vec(Array(chain[:rho]))), mean(vec(Array(chain[:alpha])))], beta, gamma)
    else
        error("Unknown model name: $model_name")
    end
end

function posterior_mean_seed(chain::Chains)
    return mean(vec(Array(chain[:seed])))
end

function write_posterior_hdf5(path::AbstractString, chain::Chains, spec::RunSpec, labels::Vector{String}, timepoints::Vector{Float64})
    samples = Array(chain)
    param_names = String.(names(chain, :parameters))
    h5open(path, "w") do h5
        h5["chains/samples"] = samples
        h5["chains/parameter_names"] = param_names
        h5["data/region_labels"] = labels
        h5["data/timepoints_train"] = timepoints
        h5["spec/model_name"] = spec.model.name
        h5["spec/transport"] = spec.model.transport
    end
    return path
end

function fit_synthetic_smoke(spec::RunSpec, L::AbstractMatrix, pathology; n_samples::Int=200)
    ignore_regions = spec.inference.ignore_seed ? spec.seeding.seed_indices : Int[]
    obs = finite_observations(
        pathology.data;
        mean_data = spec.inference.mean_data,
        ignore_regions = ignore_regions,
    )
    N = size(L, 1)
    prob = make_ode_problem(spec, Matrix{Float64}(L), ["r$(i)" for i in 1:N], pathology.timepoints)
    model = smoke_model(spec, prob, N, Float64.(pathology.timepoints), obs.values, obs)
    sampler_name = uppercase(spec.inference.sampler)
    if sampler_name == "NUTS"
        return sample(model, NUTS(spec.inference.n_warmup, spec.inference.target_acceptance; adtype = AutoReverseDiff()), n_samples; progress=false)
    elseif sampler_name == "MH"
        return sample(model, MH(), n_samples; progress=false)
    else
        error("Unsupported sampler: $(spec.inference.sampler)")
    end
end

function resolve_data_paths(spec::RunSpec, config_path::AbstractString)
    config_dir = dirname(abspath(config_path))
    network = _resolve_path(spec.data.network, config_dir)
    observations = _resolve_path(spec.data.observations, config_dir)
    return RunSpec(
        model = spec.model,
        data = DataSpec(network=network, observations=observations, region_label_style=spec.data.region_label_style),
        seeding = spec.seeding,
        inference = spec.inference,
        holdout = spec.holdout,
        run_name = spec.run_name,
    )
end

function _resolve_path(path_str::AbstractString, config_dir::AbstractString)
    if isabspath(path_str)
        return String(path_str)
    end
    candidates = [
        normpath(joinpath(config_dir, path_str)),
        normpath(joinpath(dirname(config_dir), path_str)),
        normpath(joinpath(dirname(dirname(config_dir)), path_str)),
    ]
    for candidate in candidates
        if isfile(candidate)
            return candidate
        end
    end
    return candidates[1]
end

function fit_and_save_run(spec::RunSpec; run_root::AbstractString, run_id::Union{Nothing,String}=nothing)
    transport = build_transport_operator(spec.data.network; transport=spec.model.transport)
    pathology = process_pathology(spec.data.observations; network_csv=spec.data.network)
    chain = fit_synthetic_smoke(spec, transport.L, pathology; n_samples=spec.inference.n_samples)
    params = posterior_mean_parameter_vector(chain, spec.model.name, length(transport.labels))
    seed_mean = posterior_mean_seed(chain)
    pred = simulate_trajectory(spec, transport.L, transport.labels, pathology.timepoints, params; seed_value=seed_mean)
    pred_observed = pred[1:length(transport.labels), :]

    resolved_run_id = isnothing(run_id) ? resolve_run_id(spec) : String(run_id)
    bundle_spec = portable_run_spec(spec, run_root)
    paths = initialize_run_bundle(run_root, bundle_spec; run_id=resolved_run_id, overwrite=true)
    write_posterior_hdf5(paths.posterior_path, chain, spec, transport.labels, pathology.timepoints)

    pred_df = DataFrame(region = transport.labels)
    for (j, t) in enumerate(pathology.timepoints)
        pred_df[!, string(t)] = pred_observed[:, j]
    end
    CSV.write(paths.predictions_train_path, pred_df)

    summary = posterior_summary_table(chain, parameter_names_for_model(spec.model.name, length(transport.labels)))
    CSV.write(paths.posterior_summary_path, summary)
    open(paths.diagnostics_path, "w") do io
        print(io,
            "{\"sampler\":\"$(spec.inference.sampler)\",\"n_samples\":$(spec.inference.n_samples),\"n_warmup\":$(spec.inference.n_warmup),\"status\":\"completed\"}")
    end
    return paths
end
