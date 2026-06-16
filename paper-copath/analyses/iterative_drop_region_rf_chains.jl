#!/usr/bin/env julia

using CSV
using DataFrames
using DifferentialEquations
using Distributions
using HDF5
using MCMCChains
using Statistics

const DATASETS = ["syn_app", "syn_mapt", "tau_app", "tau_mapt"]
const OBSERVATION_FILES = Dict(
    "syn_app" => "syn_pathology_app.csv",
    "syn_mapt" => "syn_pathology_mapt.csv",
    "tau_app" => "tau_pathology_app.csv",
    "tau_mapt" => "tau_pathology_mapt.csv",
)
const PARAMETERS = ["alpha", "beta", "gamma", "u0", "sigma"]
const MAIN_PARAMETERS = Set(["alpha", "beta", "gamma"])
const RHAT_CUTOFF = 1.05

function load_posterior(path::AbstractString)
    h5open(path, "r") do h5
        samples = Matrix{Float64}(read(h5["chains/samples"]))
        parameter_names = String.(read(h5["chains/parameter_names"]))
        chain_ids = Int.(read(h5["chains/chain_ids"]))
        if size(samples, 2) != length(parameter_names) && size(samples, 1) == length(parameter_names)
            samples = Matrix(permutedims(samples))
        end
        return samples, chain_ids, parameter_names
    end
end

function chain_from_retained(samples::Matrix{Float64}, chain_ids::Vector{Int}, parameter_names::Vector{String}, retained::Vector{Int})
    unique_ids = sort(retained)
    n_samples = minimum([count(==(id), chain_ids) for id in unique_ids])
    arr = Array{Float64}(undef, n_samples, size(samples, 2), length(unique_ids))
    for (j, id) in enumerate(unique_ids)
        idx = findall(==(id), chain_ids)
        arr[:, :, j] = samples[idx[1:n_samples], :]
    end
    return Chains(arr, Symbol.(parameter_names))
end

function summary_for_retained(samples::Matrix{Float64}, chain_ids::Vector{Int}, parameter_names::Vector{String}, retained::Vector{Int})
    keep = in.(chain_ids, Ref(Set(retained)))
    retained_samples = samples[keep, :]
    if length(retained) >= 2
        chain = chain_from_retained(samples, chain_ids, parameter_names, retained)
        df = DataFrame(MCMCChains.summarystats(chain))
        rename!(df, :parameters => :parameter)
        return df
    end

    rows = DataFrame(
        parameter = String[],
        mean = Float64[],
        std = Float64[],
        mcse = Float64[],
        ess_bulk = Float64[],
        ess_tail = Float64[],
        rhat = Float64[],
        ess_per_sec = Float64[],
    )
    for (j, parameter) in enumerate(parameter_names)
        values = retained_samples[:, j]
        push!(
            rows,
            (
                parameter = parameter,
                mean = mean(values),
                std = std(values),
                mcse = std(values) / sqrt(length(values)),
                ess_bulk = NaN,
                ess_tail = NaN,
                rhat = NaN,
                ess_per_sec = NaN,
            ),
        )
    end
    return rows
end

function main_max_rhat(summary::DataFrame)
    main = summary[in.(String.(summary.parameter), Ref(MAIN_PARAMETERS)), :]
    if nrow(main) == 0 || any(ismissing, main.rhat) || any(isnan, Float64.(main.rhat))
        return NaN
    end
    return maximum(Float64.(main.rhat))
end

function region_rf_rhs!(du, u, p, t)
    alpha, beta, gamma = p
    x = u[1]
    y = u[2]
    du[1] = alpha * x * (beta - y - x)
    du[2] = gamma * x
    return nothing
end

function solve_region(timepoints::Vector{Float64}, params::Vector{Float64}, u0::Float64)
    maximum(timepoints) <= 0 && return fill(u0, length(timepoints))
    prob = ODEProblem(region_rf_rhs!, [max(u0, 0.0), 0.0], (0.0, maximum(timepoints)), params)
    sol = solve(prob, Tsit5(); saveat = timepoints, abstol = 1e-8, reltol = 1e-6, maxiters = 50_000)
    if sol.retcode != ReturnCode.Success
        return nothing
    end
    pred = Array(sol)[1, :]
    all(isfinite, pred) || return nothing
    return pred
end

function chain_loglik(samples::Matrix{Float64}, chain_ids::Vector{Int}, parameter_names::Vector{String}, chain_id::Int, obs_times::Vector{Float64}, obs_values::Vector{Float64})
    idx = Dict(name => findfirst(==(name), parameter_names) for name in PARAMETERS)
    chain_samples = samples[chain_ids .== chain_id, :]
    means = Dict(name => mean(chain_samples[:, idx[name]]) for name in PARAMETERS)
    unique_times = sort(unique(obs_times))
    pred = solve_region(unique_times, [means["alpha"], means["beta"], means["gamma"]], means["u0"])
    isnothing(pred) && return -Inf
    pred_by_time = Dict(t => pred[i] for (i, t) in enumerate(unique_times))
    sigma = max(means["sigma"], 1e-12)
    return sum(logpdf.(Normal.(getindex.(Ref(pred_by_time), obs_times), sigma), obs_values))
end

function region_observations(observations::DataFrame, region::AbstractString)
    time_col = names(observations)[2]
    obs_times = Float64[]
    obs_values = Float64[]
    for row in eachrow(observations)
        value = row[Symbol(region)]
        if !ismissing(value) && isfinite(Float64(value))
            push!(obs_times, Float64(row[Symbol(time_col)]))
            push!(obs_values, Float64(value))
        end
    end
    return obs_times, obs_values
end

function original_main_max_by_region(posterior::DataFrame)
    main = posterior[in.(String.(posterior.parameter), Ref(MAIN_PARAMETERS)), :]
    return Dict(row.region_index => row.rhat_max for row in eachrow(combine(groupby(main, :region_index), :rhat => maximum => :rhat_max)))
end

function annotate_summary!(summary::DataFrame, diag)
    summary.run_id .= String(diag.run_id)
    summary.region_index .= Int(diag.region_index)
    summary.region .= String(diag.region)
    summary.rank .= Int(diag.rank)
    return summary
end

function process_dataset(project_root::AbstractString, dataset::AbstractString, out_dir::AbstractString)
    root = joinpath(project_root, "runs", "region_rf", "copath_$dataset")
    posterior = CSV.read(joinpath(root, "region_rf_posterior_summary_long.csv"), DataFrame)
    diagnostics = CSV.read(joinpath(root, "region_rf_summary.csv"), DataFrame)
    observations = CSV.read(joinpath(project_root, "paper-copath", "data", OBSERVATION_FILES[dataset]), DataFrame; missingstring = ["NA"])
    original_max = original_main_max_by_region(posterior)

    adjusted_rows = DataFrame()
    decision_rows = DataFrame()
    step_rows = DataFrame()

    for diag in eachrow(sort(diagnostics, :region_index))
        region_index = Int(diag.region_index)
        region = String(diag.region)
        run_id = String(diag.run_id)
        original_region_rows = posterior[posterior.region_index .== region_index, :]
        original_rhat = Float64(original_max[region_index])
        retained = Int[]
        dropped = Int[]
        dropped_logliks = Float64[]
        final_status = original_rhat <= RHAT_CUTOFF ? "original_converged" : "unknown"
        final_rhat = original_rhat
        final_summary = original_region_rows

        if original_rhat > RHAT_CUTOFF
            samples, chain_ids, parameter_names = load_posterior(joinpath(root, run_id, "posterior.h5"))
            retained = sort(unique(chain_ids))
            obs_times, obs_values = region_observations(observations, region)

            while true
                summary = summary_for_retained(samples, chain_ids, parameter_names, retained)
                current_rhat = main_max_rhat(summary)
                push!(
                    step_rows,
                    (
                        dataset = dataset,
                        run_id = run_id,
                        region_index = region_index,
                        region = region,
                        rank = Int(diag.rank),
                        n_retained_chains = length(retained),
                        retained_chains = join(retained, ";"),
                        main_max_rhat = current_rhat,
                    ),
                    cols = :union,
                )

                if length(retained) == 1
                    final_status = "one_chain_remaining"
                    final_rhat = current_rhat
                    final_summary = annotate_summary!(summary, diag)
                    break
                elseif isfinite(current_rhat) && current_rhat <= RHAT_CUTOFF
                    final_status = "converged_after_drop"
                    final_rhat = current_rhat
                    final_summary = annotate_summary!(summary, diag)
                    break
                end

                logliks = Dict(chain => chain_loglik(samples, chain_ids, parameter_names, chain, obs_times, obs_values) for chain in retained)
                worst_chain = sort(collect(logliks); by = x -> (x[2], x[1]))[1][1]
                push!(dropped, worst_chain)
                push!(dropped_logliks, logliks[worst_chain])
                retained = [chain for chain in retained if chain != worst_chain]
            end
        end

        append!(adjusted_rows, final_summary; cols = :union)
        push!(
            decision_rows,
            (
                dataset = dataset,
                run_id = run_id,
                region_index = region_index,
                region = region,
                rank = Int(diag.rank),
                original_main_max_rhat = original_rhat,
                final_main_max_rhat = final_rhat,
                status = final_status,
                n_dropped = length(dropped),
                dropped_chains = join(dropped, ";"),
                dropped_chain_logliks = join(round.(dropped_logliks; digits = 6), ";"),
                retained_chains = isempty(retained) ? "1;2;3;4" : join(retained, ";"),
            ),
            cols = :union,
        )
    end

    CSV.write(joinpath(out_dir, "$(dataset)_region_rf_posterior_summary_long.csv"), adjusted_rows)
    CSV.write(joinpath(out_dir, "$(dataset)_iterative_chain_drop_decisions.csv"), decision_rows)
    CSV.write(joinpath(out_dir, "$(dataset)_iterative_chain_drop_steps.csv"), step_rows)
    return adjusted_rows, decision_rows
end

function main()
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    out_dir = joinpath(project_root, "paper-copath", "results", "region_rf_iterative_drop_low_likelihood_chains")
    mkpath(out_dir)

    dataset_rows = DataFrame()
    for dataset in DATASETS
        adjusted, decisions = process_dataset(project_root, dataset, out_dir)
        main = adjusted[in.(String.(adjusted.parameter), Ref(MAIN_PARAMETERS)), :]
        region_max = combine(groupby(main, :region_index), :rhat => maximum => :max_rhat)
        append!(
            dataset_rows,
            DataFrame(
                dataset = [dataset],
                regions_originally_nonconverged = [sum(decisions.original_main_max_rhat .> RHAT_CUTOFF)],
                regions_converged_after_drop = [sum(decisions.status .== "converged_after_drop")],
                regions_one_chain_remaining = [sum(decisions.status .== "one_chain_remaining")],
                regions_originally_converged = [sum(decisions.status .== "original_converged")],
                regions_passing_final = [sum(coalesce.(region_max.max_rhat .<= RHAT_CUTOFF, false))],
                regions_total = [nrow(region_max)],
                max_rhat_final = [maximum(skipmissing(region_max.max_rhat))],
            );
            cols = :union,
        )
    end
    CSV.write(joinpath(out_dir, "iterative_chain_drop_dataset_summary.csv"), dataset_rows)
    println(out_dir)
end

main()
