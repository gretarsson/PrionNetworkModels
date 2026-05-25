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
    elseif name == "LOCAL-RF"
        return (
            alpha = truncated(Normal(0.0, 0.1), 0.0, Inf),
            beta = Normal(0.0, 1.0),
            gamma = truncated(Normal(0.0, 0.1), 0.0, Inf),
            u0 = truncated(Normal(0.0, 0.1), 0.0, Inf),
            sigma = LogNormal(0, 1),
        )
    else
        error("Unknown model name: $model_name")
    end
end

function resolve_priors(spec::RunSpec, N::Integer)
    priors = default_priors(spec.model.name, N)
    isnothing(spec.posterior_priors.source) && return priors
    return apply_posterior_priors(priors, spec.posterior_priors)
end

function apply_posterior_priors(base_priors::NamedTuple, posterior_spec::PosteriorPriorSpec)
    source = posterior_spec.source
    isnothing(source) && return base_priors

    posterior_path = posterior_source_path(source)
    posterior = load_posterior_hdf5(posterior_path)
    parameter_names = posterior.parameter_names
    prior_updates = Dict{Symbol,Any}()

    for (idx, name) in enumerate(parameter_names)
        should_update_parameter(name, posterior_spec) || continue
        key = Symbol(name)
        haskey(base_priors, key) || error("Posterior-prior parameter '$name' is not a base prior for this model")
        values = posterior.samples[:, idx]
        prior_updates[key] = posterior_prior_distribution(values, base_priors[key], posterior_spec)
    end

    requested = Set(posterior_spec.parameters)
    found = Set(name for name in parameter_names if should_update_parameter(name, posterior_spec))
    missing = setdiff(requested, found)
    isempty(missing) || error("Posterior source is missing requested parameter(s): $(join(sort(collect(missing)), ", "))")

    return merge(base_priors, (; prior_updates...))
end

function posterior_source_path(source::AbstractString)
    if isdir(source)
        path = joinpath(source, "posterior.h5")
        isfile(path) || error("Posterior-prior source directory is missing posterior.h5: $source")
        return path
    end
    isfile(source) || error("Posterior-prior source does not exist: $source")
    return String(source)
end

function should_update_parameter(name::AbstractString, spec::PosteriorPriorSpec)
    name in spec.parameters && return true
    return any(parameter_pattern_matches(name, pattern) for pattern in spec.patterns)
end

function parameter_pattern_matches(name::AbstractString, pattern::AbstractString)
    special = Set(['.', '^', '$', '+', '?', '(', ')', '{', '}', '|', '\\', '[', ']'])
    buf = IOBuffer()
    print(buf, "^")
    for ch in pattern
        if ch == '*'
            print(buf, ".*")
        elseif ch in special
            print(buf, "\\", ch)
        else
            print(buf, ch)
        end
    end
    print(buf, "\$")
    return occursin(Regex(String(take!(buf))), name)
end

function posterior_prior_distribution(values::AbstractVector, base_prior, spec::PosteriorPriorSpec)
    mu = mean(values)
    sd = max(std(values) * spec.widen, spec.min_sd)
    if prior_is_nonnegative(base_prior)
        return truncated(Normal(mu, sd), 0.0, Inf)
    end
    return Normal(mu, sd)
end

function prior_is_nonnegative(distribution; eps::Float64=1e-12)
    return !insupport(distribution, -eps) && insupport(distribution, 0.0)
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

function seed_parameter_names(spec::RunSpec)
    n_seeds = length(spec.seeding.seed_indices)
    return n_seeds == 1 ? ["seed"] : ["seed_values[$i]" for i in 1:n_seeds]
end

function parameter_names_for_model(model_name::AbstractString, N::Integer; seed_names::Vector{String}=["seed"])
    name = String(model_name)
    if name == "DIFF"
        return vcat(["rho", "sigma"], seed_names)
    elseif name == "DIFF-R"
        return vcat(["rho", "alpha"], ["beta[$i]" for i in 1:N], ["sigma"], seed_names)
    elseif name == "DIFF-RF"
        return vcat(["rho", "alpha"], ["beta[$i]" for i in 1:N], ["gamma[$i]" for i in 1:N], ["sigma"], seed_names)
    elseif name == "LOCAL-RF"
        return vcat(["alpha"], ["beta[$i]" for i in 1:N], ["gamma[$i]" for i in 1:N], ["u0[$i]" for i in 1:N], ["sigma"])
    else
        error("Unknown model name: $model_name")
    end
end

@model function inference_model(spec::RunSpec, prob::ODEProblem, N::Integer, timepoints::Vector{Float64}, obs_vals::Vector{Float64}, obs_info, priors)
    if spec.model.name == "DIFF"
        rho ~ priors.rho
        sigma ~ priors.sigma
        if length(spec.seeding.seed_indices) == 1
            seed ~ priors.seed
            seed_value = seed
        else
            seed_values ~ filldist(priors.seed, length(spec.seeding.seed_indices))
            seed_value = seed_values
        end
        params = [rho]
    elseif spec.model.name == "DIFF-R"
        rho ~ priors.rho
        alpha ~ priors.alpha
        beta ~ filldist(priors.beta, N)
        sigma ~ priors.sigma
        if length(spec.seeding.seed_indices) == 1
            seed ~ priors.seed
            seed_value = seed
        else
            seed_values ~ filldist(priors.seed, length(spec.seeding.seed_indices))
            seed_value = seed_values
        end
        params = vcat([rho, alpha], beta)
    elseif spec.model.name == "DIFF-RF"
        rho ~ priors.rho
        alpha ~ priors.alpha
        beta ~ filldist(priors.beta, N)
        gamma ~ filldist(priors.gamma, N)
        sigma ~ priors.sigma
        if length(spec.seeding.seed_indices) == 1
            seed ~ priors.seed
            seed_value = seed
        else
            seed_values ~ filldist(priors.seed, length(spec.seeding.seed_indices))
            seed_value = seed_values
        end
        params = vcat([rho, alpha], beta, gamma)
    elseif spec.model.name == "LOCAL-RF"
        alpha ~ priors.alpha
        beta ~ filldist(priors.beta, N)
        gamma ~ filldist(priors.gamma, N)
        u0 ~ filldist(priors.u0, N)
        sigma ~ priors.sigma
        params = vcat([alpha], beta, gamma)
    else
        error("Unsupported model in inference_model: $(spec.model.name)")
    end

    state0 = spec.model.name == "LOCAL-RF" ?
        initial_conditions_for_spec(spec, N; local_u0 = u0) :
        initial_conditions_for_spec(spec, N; seed_value = seed_value)

    predicted = solve(
        prob,
        Tsit5();
        u0 = state0,
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
    elseif name == "LOCAL-RF"
        beta = [mean(vec(Array(chain[Symbol("beta[$i]")]))) for i in 1:N]
        gamma = [mean(vec(Array(chain[Symbol("gamma[$i]")]))) for i in 1:N]
        return vcat([mean(vec(Array(chain[:alpha])))], beta, gamma)
    else
        error("Unknown model name: $model_name")
    end
end

function posterior_mean_seed(chain::Chains, spec::RunSpec)
    if spec.model.name == "LOCAL-RF"
        u0_names = sort(
            [String(name) for name in names(chain, :parameters) if startswith(String(name), "u0[")];
            by = name -> parse(Int, match(r"u0\[(\d+)\]", name).captures[1]),
        )
        return [mean(vec(Array(chain[Symbol(name)]))) for name in u0_names]
    end
    seed_names = seed_parameter_names(spec)
    if length(seed_names) == 1
        return mean(vec(Array(chain[Symbol(seed_names[1])])))
    end
    return [mean(vec(Array(chain[Symbol(name)]))) for name in seed_names]
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

function fit_posterior_chain(spec::RunSpec, L::AbstractMatrix, pathology; n_samples::Int=200, progress::Bool=false)
    ignore_regions = spec.inference.ignore_seed ? spec.seeding.seed_indices : Int[]
    obs = finite_observations(
        pathology.data;
        mean_data = spec.inference.mean_data,
        ignore_regions = ignore_regions,
    )
    N = size(L, 1)
    prob = make_ode_problem(spec, Matrix{Float64}(L), ["r$(i)" for i in 1:N], pathology.timepoints)
    priors = resolve_priors(spec, N)
    model = inference_model(spec, prob, N, Float64.(pathology.timepoints), obs.values, obs, priors)
    sampler_name = uppercase(spec.inference.sampler)
    if sampler_name == "NUTS"
        return sample(model, NUTS(spec.inference.n_warmup, spec.inference.target_acceptance; adtype = AutoReverseDiff()), n_samples; progress=progress)
    elseif sampler_name == "MH"
        return sample(model, MH(), n_samples; progress=progress)
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
        posterior_priors = _resolve_posterior_prior_spec(spec.posterior_priors, dirname(dirname(config_dir))),
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

function fit_and_save_run(spec::RunSpec; run_root::AbstractString, run_id::Union{Nothing,String}=nothing, progress::Bool=false)
    transport = build_transport_operator(spec.data.network; transport=spec.model.transport)
    pathology = process_pathology(spec.data.observations; network_csv=spec.data.network)
    chain = fit_posterior_chain(spec, transport.L, pathology; n_samples=spec.inference.n_samples, progress=progress)
    params = posterior_mean_parameter_vector(chain, spec.model.name, length(transport.labels))
    seed_mean = posterior_mean_seed(chain, spec)
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

    summary = posterior_summary_table(chain, parameter_names_for_model(spec.model.name, length(transport.labels); seed_names=seed_parameter_names(spec)))
    CSV.write(paths.posterior_summary_path, summary)
    open(paths.diagnostics_path, "w") do io
        print(io,
            "{\"sampler\":\"$(spec.inference.sampler)\",\"n_samples\":$(spec.inference.n_samples),\"n_warmup\":$(spec.inference.n_warmup),\"status\":\"completed\"}")
    end
    return paths
end
