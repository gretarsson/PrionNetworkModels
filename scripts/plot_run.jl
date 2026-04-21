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
    run_dir = get_arg("--run", nothing)
    isnothing(run_dir) && error("Usage: plot_run.jl --run /path/to/run_dir [--out /path/to/output_dir]")
    out_dir = get_arg("--out", nothing)
    result = plot_run_bundle(run_dir; output_dir=out_dir)
    println(result)
end

main()
