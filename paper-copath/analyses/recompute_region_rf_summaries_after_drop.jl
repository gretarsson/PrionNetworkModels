#!/usr/bin/env julia

using CSV
using DataFrames
using HDF5
using MCMCChains

const DATASETS = ["syn_app", "syn_mapt", "tau_app", "tau_mapt"]
const MAIN_PARAMETERS = Set(["alpha", "beta", "gamma"])

function posterior_chains(path::AbstractString; drop_chain::Union{Nothing,Int}=nothing)
    h5open(path, "r") do h5
        samples = Matrix{Float64}(read(h5["chains/samples"]))
        parameter_names = String.(read(h5["chains/parameter_names"]))
        chain_ids = Int.(read(h5["chains/chain_ids"]))
        if size(samples, 2) != length(parameter_names) && size(samples, 1) == length(parameter_names)
            samples = permutedims(samples)
        end
        keep = isnothing(drop_chain) ? trues(length(chain_ids)) : chain_ids .!= drop_chain
        samples = samples[keep, :]
        chain_ids = chain_ids[keep]
        unique_ids = sort(unique(chain_ids))
        n_samples = minimum([count(==(id), chain_ids) for id in unique_ids])
        arr = Array{Float64}(undef, n_samples, size(samples, 2), length(unique_ids))
        for (j, id) in enumerate(unique_ids)
            idx = findall(==(id), chain_ids)
            arr[:, :, j] = samples[idx[1:n_samples], :]
        end
        return Chains(arr, Symbol.(parameter_names))
    end
end

function summary_for_chain(chain::Chains)
    df = DataFrame(MCMCChains.summarystats(chain))
    rename!(df, :parameters => :parameter)
    return df
end

function main()
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    decision_dir = joinpath(project_root, "paper-copath", "results", "region_rf_drop_worst_chain")
    out_dir = joinpath(project_root, "paper-copath", "results", "region_rf_drop_worst_chain_mcmcchains")
    mkpath(out_dir)

    dataset_rows = DataFrame()
    for dataset in DATASETS
        root = joinpath(project_root, "runs", "region_rf", "copath_$dataset")
        original = CSV.read(joinpath(root, "region_rf_posterior_summary_long.csv"), DataFrame)
        diagnostics = CSV.read(joinpath(root, "region_rf_summary.csv"), DataFrame)
        decisions = CSV.read(joinpath(decision_dir, "$(dataset)_drop_worst_chain_decisions.csv"), DataFrame)
        decision_by_region = Dict(Int(row.region_index) => row for row in eachrow(decisions))

        rows = DataFrame()
        for diag in eachrow(sort(diagnostics, :region_index))
            region_index = Int(diag.region_index)
            decision = decision_by_region[region_index]
            if Bool(decision.was_nonconverged) && !ismissing(decision.dropped_chain)
                chain = posterior_chains(joinpath(root, String(diag.run_id), "posterior.h5"); drop_chain = Int(decision.dropped_chain))
                summary = summary_for_chain(chain)
                summary.run_id .= String(diag.run_id)
                summary.region_index .= region_index
                summary.region .= String(diag.region)
                summary.rank .= Int(diag.rank)
                append!(rows, summary; cols = :union)
            else
                append!(rows, original[original.region_index .== region_index, :]; cols = :union)
            end
        end

        CSV.write(joinpath(out_dir, "$(dataset)_region_rf_posterior_summary_long.csv"), rows)
        main_rows = rows[in.(String.(rows.parameter), Ref(MAIN_PARAMETERS)), :]
        region_max = combine(groupby(main_rows, :region_index), :rhat => maximum => :max_rhat)
        append!(
            dataset_rows,
            DataFrame(
                dataset = [dataset],
                regions_adjusted = [sum(decisions.was_nonconverged)],
                regions_passing_after = [sum(region_max.max_rhat .<= 1.05)],
                regions_total = [nrow(region_max)],
                max_rhat_after = [maximum(main_rows.rhat)],
            );
            cols = :union,
        )
    end
    CSV.write(joinpath(out_dir, "drop_worst_chain_dataset_summary.csv"), dataset_rows)
    println(out_dir)
end

main()
