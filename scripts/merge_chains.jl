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

has_flag(flag::String) = any(==(flag), ARGS)

function parse_csv_values(value::AbstractString)
    return [strip(part) for part in split(value, ",") if !isempty(strip(part))]
end

function parse_chain_ids(value::AbstractString)
    return parse.(Int, parse_csv_values(value))
end

function write_source_chain_index(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "chain,source_run_dir,archived_run_dir")
        for row in rows
            println(io, "$(row.chain),$(row.source),$(row.archive)")
        end
    end
end

function archive_source_chains!(rows; archive_root::AbstractString)
    mkpath(archive_root)
    archived_rows = map(rows) do row
        archive_path = joinpath(archive_root, basename(row.source))
        ispath(archive_path) && error("Archive destination already exists: $archive_path")
        mv(row.source, archive_path)
        (chain = row.chain, source = row.source, archive = archive_path)
    end
    return archived_rows
end

function chain_run_dir(runs_root::AbstractString, prefix::AbstractString, chain::Integer)
    run_id = "$(prefix)_C$(chain)"
    direct = joinpath(runs_root, run_id)
    isdir(direct) && return direct
    archived = joinpath(runs_root, "_source_chains", prefix, run_id)
    isdir(archived) && return archived
    return direct
end

function source_chain_rows(runs_root::AbstractString, prefix::Union{Nothing,String}, chain_labels, run_dirs)
    if isnothing(prefix)
        return [(chain = chain, source = run_dir, archive = "") for (chain, run_dir) in zip(chain_labels, run_dirs)]
    end
    return map(zip(chain_labels, run_dirs)) do (chain, run_dir)
        run_id = "$(prefix)_C$(chain)"
        source = joinpath(runs_root, run_id)
        archive = joinpath(runs_root, "_source_chains", prefix, run_id)
        archive_value = abspath(run_dir) == abspath(archive) ? archive : ""
        (chain = chain, source = source, archive = archive_value)
    end
end

function main()
    root = dirname(@__DIR__)
    runs_root = get_arg("--runs-root", joinpath(root, "runs"))
    prefix = get_arg("--prefix", nothing)
    out_run_id = get_arg("--out-run-id", prefix)
    chains_arg = get_arg("--chains", nothing)
    run_dirs_arg = get_arg("--run-dirs", nothing)
    chain_count = parse(Int, get_arg("--chain-count", "4"))
    archive_sources = has_flag("--archive-source-chains")

    if isnothing(prefix) && isnothing(run_dirs_arg)
        error("Usage: merge_chains.jl --prefix RUN_PREFIX [--chains 1,2,3 | --chain-count 4] [--out-run-id MERGED_ID] [--archive-source-chains] [--runs-root /path/to/runs]")
    end
    isnothing(out_run_id) && error("Could not determine output run id")

    if !isnothing(run_dirs_arg)
        run_dirs = parse_csv_values(run_dirs_arg)
        chain_labels = collect(1:length(run_dirs))
    elseif !isnothing(chains_arg)
        chain_labels = parse_chain_ids(chains_arg)
        run_dirs = [chain_run_dir(runs_root, prefix, i) for i in chain_labels]
    else
        chain_labels = collect(1:chain_count)
        run_dirs = [joinpath(runs_root, "$(prefix)_C$(i)") for i in chain_labels]
    end

    for run_dir in run_dirs
        isdir(run_dir) || error("Missing chain run directory: $run_dir")
    end

    paths = merge_chain_runs(run_dirs; merged_run_root = runs_root, run_id = out_run_id)
    source_rows = source_chain_rows(runs_root, prefix, chain_labels, run_dirs)
    if archive_sources
        archive_root = joinpath(runs_root, "_source_chains", out_run_id)
        source_rows = archive_source_chains!(source_rows; archive_root = archive_root)
    end
    write_source_chain_index(joinpath(paths.run_dir, "source_chains.csv"), source_rows)
    println(paths.run_dir)
end

main()
