#!/usr/bin/env julia

using CSV
using DataFrames
using PrionNetworkModels
using Random
using Statistics

function write_network(path, W, labels)
    df = DataFrame(region = labels)
    for (j, label) in enumerate(labels)
        df[!, label] = W[:, j]
    end
    CSV.write(path, df)
end

function write_observations(path, labels, timepoints, data)
    sample_ids = String[]
    times = Float64[]
    region_columns = Dict(label => Float64[] for label in labels)

    for sample in axes(data, 3), (tidx, t) in enumerate(timepoints)
        push!(sample_ids, "sample_$(sample)")
        push!(times, Float64(t))
        for (ridx, label) in enumerate(labels)
            push!(region_columns[label], Float64(data[ridx, tidx, sample]))
        end
    end

    df = DataFrame(sample_id = sample_ids, timepoint = times)
    for label in labels
        df[!, label] = region_columns[label]
    end
    CSV.write(path, df)
end

function write_observation_summary(path, labels, timepoints, data)
    summary = summarize_over_replicates(data)
    rows = DataFrame(
        region = String[],
        timepoint = Float64[],
        mean = Float64[],
        sd = Float64[],
        se = Float64[],
        n = Float64[],
    )

    for (ridx, label) in enumerate(labels), (tidx, t) in enumerate(timepoints)
        push!(rows, (
            label,
            Float64(t),
            summary.mean[ridx, tidx],
            summary.sd[ridx, tidx],
            summary.se[ridx, tidx],
            summary.n[ridx, tidx],
        ))
    end

    CSV.write(path, rows)
end

function main()
    Random.seed!(7)

    root = dirname(@__DIR__)
    labels = ["r$(i)" for i in 1:10]
    W = rand(10, 10)
    W = (W + W') ./ 2
    for i in 1:10
        W[i, i] = 0.0
    end

    network_path = joinpath(root, "data/examples/network.csv")
    observations_path = joinpath(root, "data/examples/observations.csv")
    summary_path = joinpath(root, "data/examples/observations_summary.csv")
    truth_path = joinpath(root, "data/examples/generating_parameters_diff_rf.csv")
    mkpath(dirname(network_path))

    write_network(network_path, W, labels)

    spec = load_run_spec(joinpath(root, "configs/examples/diff_rf.toml"))
    transport = build_transport_operator(network_path; transport=spec.model.transport)
    timepoints = [0.25, 0.75, 1.5, 2.5, 4.0, 6.0, 8.0, 10.0]

    rho = 0.07 + 0.03 * rand()
    alpha = 1.4 + 0.5 * rand()
    beta = 0.8 .+ 0.6 .* rand(length(labels))
    gamma = 0.08 .+ 0.08 .* rand(length(labels))
    params = vcat([rho, alpha], beta, gamma)

    full_trajectory = simulate_trajectory(spec, transport.L, transport.labels, timepoints, params; seed_value=0.5)
    x = full_trajectory[1:length(labels), :]

    n_replicates = 3
    observed3 = Array{Float64}(undef, size(x, 1), size(x, 2), n_replicates)
    for k in 1:n_replicates
        noise = 0.03 .* randn(size(x))
        observed3[:, :, k] .= max.(x .+ noise, 0.0)
    end
    write_observations(observations_path, labels, timepoints, observed3)
    write_observation_summary(summary_path, labels, timepoints, observed3)

    parameter_names = vcat(["rho", "alpha"], ["beta[$i]" for i in 1:10], ["gamma[$i]" for i in 1:10])
    summary = DataFrame(parameter = parameter_names, value = params)
    CSV.write(truth_path, summary)

    println("Created synthetic example files:")
    println(network_path)
    println(observations_path)
    println(summary_path)
    println(truth_path)
end

main()
