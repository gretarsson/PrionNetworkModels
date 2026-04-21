#!/usr/bin/env julia

using PrionNetworkModels

function get_arg(flag::String, default::Union{Nothing,String}=nothing)
    idx = findfirst(==(flag), ARGS)
    if isnothing(idx)
        return default
    end
    idx == length(ARGS) && error("Missing value for $flag")
    return ARGS[idx + 1]
end

function has_flag(flag::String)
    return flag in ARGS
end

function main()
    config_path = get_arg("--config", joinpath(@__DIR__, "..", "configs", "examples", "diff_r.toml"))
    run_id = get_arg("--run-id", nothing)
    n_samples_arg = get_arg("--samples", nothing)
    n_warmup_arg = get_arg("--warmup", nothing)

    spec = load_run_spec(config_path)
    spec = resolve_data_paths(spec, config_path)
    if n_samples_arg !== nothing || n_warmup_arg !== nothing
        inference = InferenceSpec(
            n_chains = spec.inference.n_chains,
            target_acceptance = spec.inference.target_acceptance,
            sampler = spec.inference.sampler,
            n_samples = n_samples_arg === nothing ? spec.inference.n_samples : parse(Int, n_samples_arg),
            n_warmup = n_warmup_arg === nothing ? spec.inference.n_warmup : parse(Int, n_warmup_arg),
            mean_data = spec.inference.mean_data,
            ignore_seed = spec.inference.ignore_seed,
        )
        spec = RunSpec(
            model = spec.model,
            data = spec.data,
            seeding = spec.seeding,
            inference = inference,
            holdout = spec.holdout,
            run_name = spec.run_name,
        )
    end

    paths = fit_and_save_run(spec; run_root=joinpath(dirname(@__DIR__), "runs"), run_id=run_id)
    println(paths.run_dir)
end

main()
