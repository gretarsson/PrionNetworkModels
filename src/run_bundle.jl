using Dates
using TOML
using UUIDs

Base.@kwdef struct ModelSpec
    name::String
    transport::String
    parameter_sharing::String = "independent"
end

Base.@kwdef struct DataSpec
    network::String
    observations::String
    region_label_style::String = "unspecified"
end

Base.@kwdef struct SeedingSpec
    seed_indices::Vector{Int}
    infer_seed::Bool = true
    infer_local_u0::Bool = true
    local_u0_value::Float64 = 3.364e-5
end

Base.@kwdef struct InferenceSpec
    n_chains::Int = 4
    target_acceptance::Float64 = 0.65
    sampler::String = "NUTS"
    n_samples::Int = 500
    n_warmup::Int = 500
    mean_data::Bool = false
    ignore_seed::Bool = false
end

Base.@kwdef struct HoldoutSpec
    strategy::String = "none"
    n::Int = 0
end

Base.@kwdef struct PosteriorPriorSpec
    source::Union{Nothing,String} = nothing
    parameters::Vector{String} = String[]
    patterns::Vector{String} = String[]
    widen::Float64 = 2.5
    min_sd::Float64 = 1e-6
end

Base.@kwdef struct RunSpec
    model::ModelSpec
    data::DataSpec
    seeding::SeedingSpec
    inference::InferenceSpec = InferenceSpec()
    holdout::HoldoutSpec = HoldoutSpec()
    posterior_priors::PosteriorPriorSpec = PosteriorPriorSpec()
    run_name::Union{Nothing,String} = nothing
end

Base.@kwdef struct RunBundlePaths
    run_id::String
    root_dir::String
    run_dir::String
    spec_path::String
    metadata_path::String
    posterior_path::String
    posterior_summary_path::String
    diagnostics_path::String
    predictions_train_path::String
    predictions_full_path::String
    observed_train_path::String
    observed_full_path::String
end

function load_run_spec(path::AbstractString)
    raw = TOML.parsefile(path)

    model_tbl = get(raw, "model", Dict{String,Any}())
    data_tbl = get(raw, "data", Dict{String,Any}())
    seed_tbl = get(raw, "seeding", Dict{String,Any}())
    inf_tbl = get(raw, "inference", Dict{String,Any}())
    hold_tbl = get(raw, "holdout", Dict{String,Any}())
    priors_tbl = get(raw, "priors", Dict{String,Any}())
    posterior_prior_tbl = get(priors_tbl, "posterior", Dict{String,Any}())

    model = ModelSpec(
        name = _required_string(model_tbl, "name", path),
        transport = _required_string(model_tbl, "transport", path),
        parameter_sharing = get(model_tbl, "parameter_sharing", "independent"),
    )

    data = DataSpec(
        network = _required_string(data_tbl, "network", path),
        observations = _required_string(data_tbl, "observations", path),
        region_label_style = get(data_tbl, "region_label_style", "unspecified"),
    )

    seeds_raw = get(seed_tbl, "seed_indices", nothing)
    seeds_raw === nothing && error("Missing required key [seeding].seed_indices in $path")
    seeding = SeedingSpec(
        seed_indices = Int.(seeds_raw),
        infer_seed = get(seed_tbl, "infer_seed", true),
        infer_local_u0 = get(seed_tbl, "infer_local_u0", true),
        local_u0_value = Float64(get(seed_tbl, "local_u0_value", 3.364e-5)),
    )

    inference = InferenceSpec(
        n_chains = get(inf_tbl, "n_chains", 4),
        target_acceptance = Float64(get(inf_tbl, "target_acceptance", 0.65)),
        sampler = get(inf_tbl, "sampler", "NUTS"),
        n_samples = get(inf_tbl, "n_samples", 500),
        n_warmup = get(inf_tbl, "n_warmup", 500),
        mean_data = get(inf_tbl, "mean_data", false),
        ignore_seed = get(inf_tbl, "ignore_seed", false),
    )

    holdout = HoldoutSpec(
        strategy = get(hold_tbl, "strategy", "none"),
        n = get(hold_tbl, "n", 0),
    )

    posterior_priors = PosteriorPriorSpec(
        source = get(posterior_prior_tbl, "source", nothing),
        parameters = String.(get(posterior_prior_tbl, "parameters", String[])),
        patterns = String.(get(posterior_prior_tbl, "patterns", String[])),
        widen = Float64(get(posterior_prior_tbl, "widen", 2.5)),
        min_sd = Float64(get(posterior_prior_tbl, "min_sd", 1e-6)),
    )

    return RunSpec(
        model = model,
        data = data,
        seeding = seeding,
        inference = inference,
        holdout = holdout,
        posterior_priors = posterior_priors,
        run_name = get(raw, "run_name", nothing),
    )
end

function spec_to_dict(spec::RunSpec)
    out = Dict{String, Any}(
        "model" => Dict(
            "name" => spec.model.name,
            "transport" => spec.model.transport,
            "parameter_sharing" => spec.model.parameter_sharing,
        ),
        "data" => Dict(
            "network" => spec.data.network,
            "observations" => spec.data.observations,
            "region_label_style" => spec.data.region_label_style,
        ),
        "seeding" => Dict(
            "seed_indices" => spec.seeding.seed_indices,
            "infer_seed" => spec.seeding.infer_seed,
            "infer_local_u0" => spec.seeding.infer_local_u0,
            "local_u0_value" => spec.seeding.local_u0_value,
        ),
        "inference" => Dict(
            "n_chains" => spec.inference.n_chains,
            "target_acceptance" => spec.inference.target_acceptance,
            "sampler" => spec.inference.sampler,
            "n_samples" => spec.inference.n_samples,
            "n_warmup" => spec.inference.n_warmup,
            "mean_data" => spec.inference.mean_data,
            "ignore_seed" => spec.inference.ignore_seed,
        ),
        "holdout" => Dict(
            "strategy" => spec.holdout.strategy,
            "n" => spec.holdout.n,
        ),
    )
    if !isnothing(spec.posterior_priors.source)
        out["priors"] = Dict(
            "posterior" => Dict(
                "source" => spec.posterior_priors.source,
                "parameters" => spec.posterior_priors.parameters,
                "patterns" => spec.posterior_priors.patterns,
                "widen" => spec.posterior_priors.widen,
                "min_sd" => spec.posterior_priors.min_sd,
            ),
        )
    end
    if !isnothing(spec.run_name)
        out["run_name"] = spec.run_name
    end
    return out
end

function portable_run_spec(spec::RunSpec, run_root::AbstractString)
    project_root = dirname(abspath(run_root))
    return RunSpec(
        model = spec.model,
        data = DataSpec(
            network = _portable_path(spec.data.network, project_root),
            observations = _portable_path(spec.data.observations, project_root),
            region_label_style = spec.data.region_label_style,
        ),
        seeding = spec.seeding,
        inference = spec.inference,
        holdout = spec.holdout,
        posterior_priors = _portable_posterior_prior_spec(spec.posterior_priors, project_root),
        run_name = spec.run_name,
    )
end

function resolve_bundle_spec_paths(spec::RunSpec, run_dir::AbstractString)
    project_root = _find_project_root(run_dir)
    return RunSpec(
        model = spec.model,
        data = DataSpec(
            network = _resolve_bundle_path(spec.data.network, project_root),
            observations = _resolve_bundle_path(spec.data.observations, project_root),
            region_label_style = spec.data.region_label_style,
        ),
        seeding = spec.seeding,
        inference = spec.inference,
        holdout = spec.holdout,
        posterior_priors = _resolve_posterior_prior_spec(spec.posterior_priors, project_root),
        run_name = spec.run_name,
    )
end

function _find_project_root(path::AbstractString)
    current = abspath(path)
    if isfile(current)
        current = dirname(current)
    end
    while true
        isfile(joinpath(current, "Project.toml")) && return current
        parent = dirname(current)
        parent == current && break
        current = parent
    end
    return dirname(dirname(abspath(path)))
end

function _portable_posterior_prior_spec(spec::PosteriorPriorSpec, project_root::AbstractString)
    isnothing(spec.source) && return spec
    return PosteriorPriorSpec(
        source = _portable_path(spec.source, project_root),
        parameters = spec.parameters,
        patterns = spec.patterns,
        widen = spec.widen,
        min_sd = spec.min_sd,
    )
end

function _resolve_posterior_prior_spec(spec::PosteriorPriorSpec, project_root::AbstractString)
    isnothing(spec.source) && return spec
    source = isabspath(spec.source) ? spec.source : normpath(joinpath(project_root, spec.source))
    return PosteriorPriorSpec(
        source = source,
        parameters = spec.parameters,
        patterns = spec.patterns,
        widen = spec.widen,
        min_sd = spec.min_sd,
    )
end

function resolve_run_id(spec::RunSpec; timestamp::DateTime=now())
    prefix = isnothing(spec.run_name) ? Dates.format(timestamp, dateformat"yyyy-mm-dd_HHMMSS") : slugify(spec.run_name)
    model = slugify(spec.model.name)
    transport = slugify(spec.model.transport)
    sharing = slugify(spec.model.parameter_sharing)
    return string(prefix, "_", model, "_", transport, "_", sharing)
end

function bundle_paths(root_dir::AbstractString, run_id::AbstractString)
    run_dir = joinpath(root_dir, run_id)
    return RunBundlePaths(
        run_id = run_id,
        root_dir = String(root_dir),
        run_dir = run_dir,
        spec_path = joinpath(run_dir, "spec.toml"),
        metadata_path = joinpath(run_dir, "metadata.json"),
        posterior_path = joinpath(run_dir, "posterior.h5"),
        posterior_summary_path = joinpath(run_dir, "posterior_summary.csv"),
        diagnostics_path = joinpath(run_dir, "diagnostics.json"),
        predictions_train_path = joinpath(run_dir, "predictions_train.csv"),
        predictions_full_path = joinpath(run_dir, "predictions_full.csv"),
        observed_train_path = joinpath(run_dir, "observed_train.csv"),
        observed_full_path = joinpath(run_dir, "observed_full.csv"),
    )
end

function initialize_run_bundle(root_dir::AbstractString, spec::RunSpec; run_id::Union{Nothing,String}=nothing, overwrite::Bool=false)
    resolved_run_id = isnothing(run_id) ? resolve_run_id(spec) : String(run_id)
    paths = bundle_paths(root_dir, resolved_run_id)

    if isdir(paths.run_dir) && !overwrite
        error("Run directory already exists: $(paths.run_dir)")
    end

    if isdir(paths.run_dir) && overwrite
        rm(paths.run_dir; force = true, recursive = true)
    end

    mkpath(paths.run_dir)
    _write_toml(paths.spec_path, spec_to_dict(spec))

    metadata = Dict(
        "run_id" => paths.run_id,
        "bundle_version" => 1,
        "created_at" => string(now()),
        "uuid" => string(uuid4()),
        "files" => Dict(
            "spec" => basename(paths.spec_path),
            "metadata" => basename(paths.metadata_path),
            "posterior" => basename(paths.posterior_path),
            "posterior_summary" => basename(paths.posterior_summary_path),
            "diagnostics" => basename(paths.diagnostics_path),
            "predictions_train" => basename(paths.predictions_train_path),
            "predictions_full" => basename(paths.predictions_full_path),
        ),
        "spec" => spec_to_dict(spec),
    )
    _write_json(paths.metadata_path, metadata)

    return paths
end

function slugify(s::AbstractString)
    buf = IOBuffer()
    prev_dash = false
    for ch in lowercase(s)
        if isletter(ch) || isnumeric(ch)
            print(buf, ch)
            prev_dash = false
        elseif !prev_dash
            print(buf, '-')
            prev_dash = true
        end
    end
    return strip(String(take!(buf)), '-')
end

function _required_string(tbl::AbstractDict, key::AbstractString, path::AbstractString)
    haskey(tbl, key) || error("Missing required key [$key] in $path")
    return String(tbl[key])
end

function _portable_path(path_str::AbstractString, project_root::AbstractString)
    absolute = abspath(path_str)
    root = normpath(project_root)
    if _is_subpath(absolute, root)
        return relpath(absolute, root)
    end
    return String(path_str)
end

function _resolve_bundle_path(path_str::AbstractString, project_root::AbstractString)
    if isabspath(path_str)
        return String(path_str)
    end
    resolved = normpath(joinpath(project_root, path_str))
    if isfile(resolved)
        return resolved
    end

    legacy_prefix = "data/paper/"
    if startswith(path_str, legacy_prefix)
        migrated = normpath(joinpath(project_root, "paper-rf", "data", path_str[length(legacy_prefix)+1:end]))
        isfile(migrated) && return migrated
    end

    return resolved
end

function _is_subpath(path::AbstractString, root::AbstractString)
    norm_path = normpath(path)
    norm_root = normpath(root)
    return norm_path == norm_root || startswith(norm_path, norm_root * Base.Filesystem.path_separator)
end

function _write_toml(path::AbstractString, obj::AbstractDict)
    open(path, "w") do io
        TOML.print(io, obj)
    end
    return path
end

function _write_json(path::AbstractString, obj)
    open(path, "w") do io
        print(io, _json(obj))
    end
    return path
end

function _json(x::Nothing)
    return "null"
end

function _json(x::Bool)
    return x ? "true" : "false"
end

function _json(x::Integer)
    return string(x)
end

function _json(x::AbstractFloat)
    return isfinite(x) ? string(x) : error("JSON writer does not support non-finite float values")
end

function _json(x::AbstractString)
    escaped = replace(x,
        "\\" => "\\\\",
        "\"" => "\\\"",
        "\n" => "\\n",
        "\r" => "\\r",
        "\t" => "\\t",
    )
    return "\"" * escaped * "\""
end

function _json(xs::AbstractVector)
    return "[" * join((_json(x) for x in xs), ",") * "]"
end

function _json(d::AbstractDict)
    parts = String[]
    for key in sort!(String.(collect(keys(d))))
        push!(parts, _json(key) * ":" * _json(d[key]))
    end
    return "{" * join(parts, ",") * "}"
end

function _json(x)
    error("Unsupported JSON value of type $(typeof(x))")
end
