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

function main()
    root = dirname(@__DIR__)
    runs_root = get_arg("--runs-root", joinpath(root, "runs"))
    prefix = get_arg("--prefix", nothing)
    out_run_id = get_arg("--out-run-id", prefix)
    chain_count = parse(Int, get_arg("--chain-count", "4"))

    isnothing(prefix) && error("Usage: merge_chains.jl --prefix RUN_PREFIX [--out-run-id MERGED_ID] [--chain-count 4] [--runs-root /path/to/runs]")
    isnothing(out_run_id) && error("Could not determine output run id")

    run_dirs = [joinpath(runs_root, "$(prefix)_C$(i)") for i in 1:chain_count]
    for run_dir in run_dirs
        isdir(run_dir) || error("Missing chain run directory: $run_dir")
    end

    paths = merge_chain_runs(run_dirs; merged_run_root = runs_root, run_id = out_run_id)
    println(paths.run_dir)
end

main()
