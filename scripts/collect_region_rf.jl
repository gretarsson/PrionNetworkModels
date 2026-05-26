#!/usr/bin/env julia

using CSV
using DataFrames

function get_arg(flag::String, default::Union{Nothing,String}=nothing)
    idx = findfirst(==(flag), ARGS)
    isnothing(idx) && return default
    idx == length(ARGS) && error("Missing value for $flag")
    return ARGS[idx + 1]
end

function diagnostics_paths(root::AbstractString)
    paths = String[]
    for (dirpath, _, filenames) in walkdir(root)
        "diagnostics.csv" in filenames && push!(paths, joinpath(dirpath, "diagnostics.csv"))
    end
    sort!(paths)
    return paths
end

function posterior_summary_paths(root::AbstractString)
    paths = String[]
    for (dirpath, _, filenames) in walkdir(root)
        if "diagnostics.csv" in filenames && "posterior_summary.csv" in filenames
            push!(paths, joinpath(dirpath, "posterior_summary.csv"))
        end
    end
    sort!(paths)
    return paths
end

function main()
    root = get_arg("--root", nothing)
    isnothing(root) && error("Usage: collect_region_rf.jl --root runs/region_rf/DATASET [--out summary.csv]")
    out = get_arg("--out", joinpath(root, "region_rf_summary.csv"))
    posterior_out = get_arg("--posterior-out", joinpath(root, "region_rf_posterior_summary_long.csv"))

    diag_files = diagnostics_paths(root)
    isempty(diag_files) && error("No diagnostics.csv files found under $root")
    diagnostics_tables = CSV.read.(diag_files, Ref(DataFrame))
    diagnostics = reduce((a, b) -> vcat(a, b; cols = :union), diagnostics_tables)
    sort!(diagnostics, [:region_index])
    CSV.write(out, diagnostics)

    posterior_rows = DataFrame()
    for path in posterior_summary_paths(root)
        diag_path = joinpath(dirname(path), "diagnostics.csv")
        diag = CSV.read(diag_path, DataFrame)
        posterior = CSV.read(path, DataFrame)
        posterior[!, :run_id] = fill(diag.run_id[1], nrow(posterior))
        posterior[!, :region_index] = fill(diag.region_index[1], nrow(posterior))
        posterior[!, :region] = fill(diag.region[1], nrow(posterior))
        posterior[!, :rank] = fill(diag.rank[1], nrow(posterior))
        append!(posterior_rows, posterior; cols = :union)
    end
    if nrow(posterior_rows) > 0
        sort!(posterior_rows, [:region_index, :parameter])
        CSV.write(posterior_out, posterior_rows)
    end

    println(out)
    println(posterior_out)
end

main()
