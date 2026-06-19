#!/usr/bin/env julia

using CSV
using DataFrames
using Dates
using DifferentialEquations
using HDF5
using PathoSpread
using Serialization
using Statistics
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SYN_ROOT = get(ENV, "SYNUCLEIN_SPREAD_ROOT", "/Users/gretarsson/Desktop/synuclein_spread")
const RUN_ROOT = joinpath(PROJECT_ROOT, "runs")
const WAIC_SAMPLES = parse(Int, get(ENV, "WAIC_POSTERIOR_SAMPLES", "300"))

const IMPORTS = [
    (
        source = "simulations/u0_DIFF_RETRO.jls",
        run_id = "striatum_DIFF_RETRO",
        model = "DIFF",
        transport = "retrograde",
        holdout = "none",
    ),
    (
        source = "simulations/u0_DIFF_EUCL.jls",
        run_id = "striatum_DIFF_EUCL",
        model = "DIFF",
        transport = "euclidean",
        holdout = "none",
    ),
    (
        source = "simulations/u0_DIFF_ANTERO.jls",
        run_id = "striatum_DIFF_ANTERO",
        model = "DIFF",
        transport = "anterograde",
        holdout = "none",
    ),
    (
        source = "simulations/u0_DIFF_BIDIR.jls",
        run_id = "striatum_DIFF_BIDIR",
        model = "DIFF",
        transport = "bidirectional",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFG_RETRO.jls",
        run_id = "striatum_DIFF-R_RETRO",
        model = "DIFF-R",
        transport = "retrograde",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFG_EUCL.jls",
        run_id = "striatum_DIFF-R_EUCL",
        model = "DIFF-R",
        transport = "euclidean",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFG_ANTERO.jls",
        run_id = "striatum_DIFF-R_ANTERO",
        model = "DIFF-R",
        transport = "anterograde",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFG_BIDIR.jls",
        run_id = "striatum_DIFF-R_BIDIR",
        model = "DIFF-R",
        transport = "bidirectional",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFGA_RETRO.jls",
        run_id = "striatum_DIFF-RF_RETRO",
        model = "DIFF-RF",
        transport = "retrograde",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFGA_EUCL.jls",
        run_id = "striatum_DIFF-RF_EUCL",
        model = "DIFF-RF",
        transport = "euclidean",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFGA_ANTERO.jls",
        run_id = "striatum_DIFF-RF_ANTERO",
        model = "DIFF-RF",
        transport = "anterograde",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFGA_BIDIR.jls",
        run_id = "striatum_DIFF-RF_BIDIR",
        model = "DIFF-RF",
        transport = "bidirectional",
        holdout = "none",
    ),
    (
        source = "simulations/DIFFG_T1.jls",
        run_id = "striatum_DIFF-R_RETRO_T-1",
        model = "DIFF-R",
        transport = "retrograde",
        holdout = "leave_final_timepoint_out",
    ),
    (
        source = "simulations/DIFFGA_T1.jls",
        run_id = "striatum_DIFF-RF_RETRO_T-1",
        model = "DIFF-RF",
        transport = "retrograde",
        holdout = "leave_final_timepoint_out",
    ),
]

function normalized_parameter_names(priors)
    names = String.(collect(keys(priors)))
    return [name == "σ" ? "sigma" : name for name in names]
end

function flattened_parameter_samples(chain, n_parameters::Integer)
    value_data = chain.value.data[:, 1:n_parameters, :]
    n_iter, _, n_chains = size(value_data)
    samples = Matrix{Float64}(undef, n_iter * n_chains, n_parameters)
    chain_ids = Vector{Int}(undef, n_iter * n_chains)
    cursor = 1
    for chain_idx in 1:n_chains
        range = cursor:(cursor + n_iter - 1)
        samples[range, :] = value_data[:, :, chain_idx]
        chain_ids[range] .= chain_idx
        cursor += n_iter
    end
    return samples, chain_ids
end

function write_spec(path::AbstractString, item, inf, n_samples::Integer, n_chains::Integer)
    seed_idx = get(inf, "seed_idx", [74])
    if seed_idx isa Integer
        seed_idx = [Int(seed_idx)]
    end
    spec = Dict{String,Any}(
        "model" => Dict(
            "name" => item.model,
            "transport" => item.transport,
            "parameter_sharing" => "independent",
        ),
        "data" => Dict(
            "observations" => "paper-rf/data/striatum/observations.csv",
            "network" => "paper-rf/data/striatum/network.csv",
            "region_label_style" => "unspecified",
        ),
        "seeding" => Dict(
            "seed_indices" => Int.(seed_idx),
            "infer_seed" => true,
            "infer_local_u0" => true,
            "local_u0_value" => 3.364e-5,
        ),
        "inference" => Dict(
            "n_samples" => n_samples,
            "n_warmup" => 1000,
            "mean_data" => false,
            "ignore_seed" => false,
            "n_chains" => n_chains,
            "target_acceptance" => 0.65,
            "sampler" => "NUTS",
        ),
        "holdout" => Dict(
            "strategy" => item.holdout,
            "n" => item.holdout == "none" ? 0 : 1,
        ),
    )
    open(path, "w") do io
        TOML.print(io, spec)
    end
end

function write_summary(path::AbstractString, samples::AbstractMatrix, parameter_names::Vector{String})
    rows = DataFrame(
        parameter = parameter_names,
        mean = [mean(samples[:, i]) for i in axes(samples, 2)],
        sd = [std(samples[:, i]) for i in axes(samples, 2)],
    )
    CSV.write(path, rows)
end

function posterior_mode_prediction(inf, timepoints::AbstractVector)
    chain = inf["chain"]
    prior_names = collect(keys(inf["priors"]))
    n_parameters = findfirst(==("σ"), prior_names) - 1
    n_chain_parameters = length(chain.info[1])
    _, mode_index = findmax(chain[:lp])
    mode_pars = Array(chain[mode_index[1], 1:n_chain_parameters, mode_index[2]])
    p = mode_pars[1:n_parameters]

    u0 = copy(inf["u0"])
    seed = inf["seed_idx"]
    if get(inf, "bayesian_seed", false)
        chain_parameter_names = names(chain, :parameters)
        seed_chain_indices = sort(findall(name -> startswith(String(name), "seed"), chain_parameter_names))
        isempty(seed_chain_indices) && error("No seed parameters found in chain for $(inf["ode"])")
        if seed isa Int
            u0[seed] = chain.value[mode_index[1], seed_chain_indices[1], mode_index[2]]
        else
            for (i, seed_idx) in enumerate(seed)
                u0[seed_idx] = chain.value[mode_index[1], seed_chain_indices[i], mode_index[2]]
            end
        end
    end

    ode = PathoSpread.odes[inf["ode"]]
    factors = [1.0 for _ in 1:n_parameters]
    prob = PathoSpread.make_ode_problem(
        ode;
        labels = inf["labels"],
        Ltuple = inf["L"],
        factors = factors,
        u0 = u0,
        timepoints = timepoints,
        seed_indices = inf["seed_idx"],
    )
    sol = solve(prob, Tsit5(); p = p, u0 = u0, saveat = timepoints, abstol = 1e-9, reltol = 1e-6)
    return Array(sol[inf["sol_idxs"], :]), mode_index
end

function write_prediction_table(path::AbstractString, labels, timepoints, pred)
    pred_df = DataFrame(region = String.(labels))
    for (j, timepoint) in enumerate(timepoints)
        pred_df[!, string(Float64(timepoint))] = pred[:, j]
    end
    CSV.write(path, pred_df)
end

function write_mode_predictions(path::AbstractString, inf)
    pred, mode_index = posterior_mode_prediction(inf, inf["timepoints"])
    write_prediction_table(path, inf["labels"], inf["timepoints"], pred)
    return mode_index
end

function write_holdout_mode_predictions(path::AbstractString, inf)
    _, full_timepoints = PathoSpread.process_pathology(
        joinpath(SYN_ROOT, "data", "total_path.csv");
        W_csv = joinpath(SYN_ROOT, "data", "W_labeled_filtered.csv"),
    )
    pred, mode_index = posterior_mode_prediction(inf, full_timepoints)
    write_prediction_table(path, inf["labels"], full_timepoints, pred)
    return mode_index
end

function write_dense_mode_predictions(path::AbstractString, inf; max_timepoint::Union{Nothing,Real}=nothing)
    max_t = isnothing(max_timepoint) ? maximum(Float64.(inf["timepoints"])) : Float64(max_timepoint)
    dense_timepoints = collect(range(0.0, max_t; length=300))
    pred, mode_index = posterior_mode_prediction(inf, dense_timepoints)
    write_prediction_table(path, inf["labels"], dense_timepoints, pred)
    return mode_index
end

function write_waic_tables(run_dir::AbstractString, inf)
    waic, se_waic, waic_i, lppd, p_waic, n_used = PathoSpread.compute_waic(
        inf;
        S = WAIC_SAMPLES,
        group_cells = false,
        ignore_seed = true,
    )
    CSV.write(joinpath(run_dir, "waic_summary.csv"), DataFrame(
        waic = [waic],
        se_waic = [se_waic],
        lppd = [lppd],
        p_waic = [p_waic],
        n_used = [n_used],
        posterior_samples = [WAIC_SAMPLES],
        group_cells = [false],
        ignore_seed = [true],
    ))
    CSV.write(joinpath(run_dir, "waic_pointwise.csv"), DataFrame(index = eachindex(waic_i), waic_i = waic_i))
end

function read_cached_waic_values(prefix::AbstractString; exclude::Union{Nothing,Regex}=nothing)
    cache_dir = joinpath(SYN_ROOT, "results", "waic_cache")
    files = sort(filter(f -> occursin(Regex("^$(prefix)_\\d+_waic\\.jls\$"), basename(f)), readdir(cache_dir; join=true)))
    rows = DataFrame(source = String[], waic = Float64[])
    for file in files
        base = replace(basename(file), "_waic.jls" => "")
        !isnothing(exclude) && occursin(exclude, base) && continue
        value = Serialization.deserialize(file)
        isfinite(value) || continue
        value < 0 || continue
        push!(rows, (base, value))
    end
    return rows
end

function write_null_waic_bundles()
    mkpath(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_connectivity_nulls"))
    mkpath(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_seed_nulls"))
    connectivity = read_cached_waic_values("DIFFGA_shuffle")
    seeding = read_cached_waic_values("DIFFGA_seed"; exclude = r"DIFFGA_seed_74")
    CSV.write(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_connectivity_nulls", "waic_values.csv"), connectivity)
    CSV.write(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_seed_nulls", "waic_values.csv"), seeding)
    write_json_object(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_connectivity_nulls", "diagnostics.json"), Dict(
        "status" => "imported_waic_cache_from_synuclein_spread",
        "source" => joinpath(SYN_ROOT, "results", "waic_cache", "DIFFGA_shuffle_*_waic.jls"),
        "n_values" => nrow(connectivity),
    ))
    write_json_object(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_seed_nulls", "diagnostics.json"), Dict(
        "status" => "imported_waic_cache_from_synuclein_spread",
        "source" => joinpath(SYN_ROOT, "results", "waic_cache", "DIFFGA_seed_*_waic.jls"),
        "n_values" => nrow(seeding),
    ))
end

function escape_json_string(value::AbstractString)
    return replace(value, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
end

function json_value(value)
    if value isa AbstractString
        return "\"" * escape_json_string(value) * "\""
    elseif value isa AbstractVector
        return "[" * join(json_value.(value), ", ") * "]"
    elseif value isa Number || value isa Bool
        return string(value)
    else
        return "\"" * escape_json_string(string(value)) * "\""
    end
end

function write_json_object(path::AbstractString, dict::Dict{String,<:Any})
    keys_sorted = sort(collect(keys(dict)))
    open(path, "w") do io
        println(io, "{")
        for (i, key) in enumerate(keys_sorted)
            suffix = i == length(keys_sorted) ? "" : ","
            println(io, "  \"$(escape_json_string(key))\": $(json_value(dict[key]))$suffix")
        end
        println(io, "}")
    end
end

function import_one(item)
    source_path = joinpath(SYN_ROOT, item.source)
    isfile(source_path) || error("Missing source inference: $source_path")
    run_dir = joinpath(RUN_ROOT, item.run_id)
    mkpath(run_dir)
    for stale in ("metadata.json", "source_chains.csv")
        stale_path = joinpath(run_dir, stale)
        isfile(stale_path) && rm(stale_path)
    end

    inf = PathoSpread.load_inference(source_path)
    parameter_names = normalized_parameter_names(inf["priors"])
    samples, chain_ids = flattened_parameter_samples(inf["chain"], length(parameter_names))
    n_samples_per_chain = size(inf["chain"], 1)
    n_chains = size(inf["chain"], 3)

    h5open(joinpath(run_dir, "posterior.h5"), "w") do h5
        h5["chains/samples"] = samples
        h5["chains/parameter_names"] = parameter_names
        h5["chains/chain_ids"] = chain_ids
        h5["data/region_labels"] = String.(inf["labels"])
        h5["data/timepoints_train"] = Float64.(inf["timepoints"])
        h5["spec/model_name"] = item.model
        h5["spec/transport"] = item.transport
    end

    write_spec(joinpath(run_dir, "spec.toml"), item, inf, n_samples_per_chain, n_chains)
    write_summary(joinpath(run_dir, "posterior_summary.csv"), samples, parameter_names)
    mode_index = item.holdout == "none" ?
        write_mode_predictions(joinpath(run_dir, "predictions_train.csv"), inf) :
        write_holdout_mode_predictions(joinpath(run_dir, "predictions_train.csv"), inf)
    dense_mode_index = item.holdout == "none" ?
        write_dense_mode_predictions(joinpath(run_dir, "predictions_mode_dense.csv"), inf) :
        write_dense_mode_predictions(joinpath(run_dir, "predictions_mode_dense.csv"), inf; max_timepoint=9.0)
    write_waic_tables(run_dir, inf)
    write_json_object(joinpath(run_dir, "metadata.json"), Dict(
        "bundle_version" => 1,
        "created_at" => string(now()),
        "run_id" => item.run_id,
        "source_format" => "synuclein_spread_jls",
        "source" => source_path,
        "model" => item.model,
        "transport" => item.transport,
        "prediction_summary" => item.holdout == "none" ? "posterior_mode_train_timepoints_and_dense_grid" : "posterior_mode_full_timepoints_and_dense_grid",
    ))
    write_json_object(joinpath(run_dir, "diagnostics.json"), Dict(
        "status" => "imported_from_synuclein_spread",
        "source" => source_path,
        "imported_at" => string(now()),
        "n_samples_total" => size(samples, 1),
        "n_chains" => n_chains,
        "prediction_summary" => item.holdout == "none" ? "posterior_mode_train_timepoints_and_dense_grid" : "posterior_mode_full_timepoints_and_dense_grid",
        "posterior_mode_chain" => isnothing(mode_index) ? "" : mode_index[2],
        "posterior_mode_iteration" => isnothing(mode_index) ? "" : mode_index[1],
        "posterior_dense_mode_chain" => isnothing(dense_mode_index) ? "" : dense_mode_index[2],
        "posterior_dense_mode_iteration" => isnothing(dense_mode_index) ? "" : dense_mode_index[1],
    ))
    write_json_object(joinpath(run_dir, "source_synuclein_spread.json"), Dict(
        "source" => source_path,
        "source_ode" => string(inf["ode"]),
        "source_timepoints" => Float64.(inf["timepoints"]),
        "run_id" => item.run_id,
    ))
    println("Imported $(item.source) -> runs/$(item.run_id)")
end

function main()
    mkpath(RUN_ROOT)
    for item in IMPORTS
        import_one(item)
    end
    write_null_waic_bundles()
end

main()
