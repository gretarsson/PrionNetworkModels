function posterior_summary_table(samples::AbstractMatrix, parameter_names::Vector{String})
    out_names = String[]
    means = Float64[]
    sds = Float64[]
    for (j, name) in enumerate(parameter_names)
        values = samples[:, j]
        push!(out_names, name)
        push!(means, mean(values))
        push!(sds, std(values))
    end
    return DataFrame(parameter = out_names, mean = means, sd = sds)
end

function posterior_mean_parameter_vector(samples::AbstractMatrix, parameter_names::Vector{String}, model_name::AbstractString, N::Integer)
    name_to_idx = Dict(name => idx for (idx, name) in enumerate(parameter_names))
    if model_name == "DIFF"
        return [mean(samples[:, name_to_idx["rho"]])]
    elseif model_name == "DIFF-R"
        beta = [mean(samples[:, name_to_idx["beta[$i]"]]) for i in 1:N]
        return vcat([mean(samples[:, name_to_idx["rho"]]), mean(samples[:, name_to_idx["alpha"]])], beta)
    elseif model_name == "DIFF-RF"
        beta = [mean(samples[:, name_to_idx["beta[$i]"]]) for i in 1:N]
        gamma = [mean(samples[:, name_to_idx["gamma[$i]"]]) for i in 1:N]
        return vcat([mean(samples[:, name_to_idx["rho"]]), mean(samples[:, name_to_idx["alpha"]])], beta, gamma)
    else
        error("Unknown model name: $model_name")
    end
end

function posterior_mean_seed(samples::AbstractMatrix, parameter_names::Vector{String})
    idx = findfirst(==("seed"), parameter_names)
    if !isnothing(idx)
        return mean(samples[:, idx])
    end

    seed_idxs = findall(name -> startswith(name, "seed_values["), parameter_names)
    isempty(seed_idxs) && error("Posterior samples do not contain seed parameters")
    sort!(seed_idxs; by = idx -> parameter_names[idx])
    return [mean(samples[:, idx]) for idx in seed_idxs]
end

function load_posterior_hdf5(path::AbstractString)
    h5open(path, "r") do h5
        samples = Matrix{Float64}(read(h5["chains/samples"]))
        parameter_names = String.(read(h5["chains/parameter_names"]))
        chain_ids = haskey(h5, "chains/chain_ids") ? Int.(read(h5["chains/chain_ids"])) : collect(ones(Int, size(samples, 1)))
        return (samples = samples, parameter_names = parameter_names, chain_ids = chain_ids)
    end
end

function write_posterior_hdf5(path::AbstractString, samples::AbstractMatrix, parameter_names::Vector{String}, spec::RunSpec, labels::Vector{String}, timepoints::Vector{Float64}; chain_ids::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    h5open(path, "w") do h5
        h5["chains/samples"] = Matrix{Float64}(samples)
        h5["chains/parameter_names"] = parameter_names
        if !isnothing(chain_ids)
            h5["chains/chain_ids"] = Int.(collect(chain_ids))
        end
        h5["data/region_labels"] = labels
        h5["data/timepoints_train"] = timepoints
        h5["spec/model_name"] = spec.model.name
        h5["spec/transport"] = spec.model.transport
    end
    return path
end

function merge_chain_runs(run_dirs::Vector{String}; merged_run_root::AbstractString, run_id::AbstractString)
    isempty(run_dirs) && error("No run directories were provided for merging")

    raw_specs = [load_run_spec(joinpath(run_dir, "spec.toml")) for run_dir in run_dirs]
    specs = [resolve_bundle_spec_paths(spec, run_dir) for (spec, run_dir) in zip(raw_specs, run_dirs)]
    first_spec = first(specs)
    first_spec_dict = spec_to_dict(first_spec)
    for spec in specs[2:end]
        spec_to_dict(spec) == first_spec_dict || error("Run specs do not match and cannot be merged")
    end

    posterior_blocks = [load_posterior_hdf5(joinpath(run_dir, "posterior.h5")) for run_dir in run_dirs]
    parameter_names = posterior_blocks[1].parameter_names
    for block in posterior_blocks[2:end]
        block.parameter_names == parameter_names || error("Posterior parameter names differ across chain runs")
    end

    samples = reduce(vcat, [block.samples for block in posterior_blocks])
    chain_ids = reduce(vcat, [fill(i, size(block.samples, 1)) for (i, block) in enumerate(posterior_blocks)])

    transport = build_transport_operator(first_spec.data.network; transport = first_spec.model.transport)
    pathology = process_pathology(first_spec.data.observations; network_csv = first_spec.data.network)
    params = posterior_mean_parameter_vector(samples, parameter_names, first_spec.model.name, length(transport.labels))
    seed_mean = posterior_mean_seed(samples, parameter_names)
    pred = simulate_trajectory(first_spec, transport.L, transport.labels, pathology.timepoints, params; seed_value = seed_mean)
    pred_observed = pred[1:length(transport.labels), :]

    bundle_spec = portable_run_spec(first_spec, merged_run_root)
    paths = initialize_run_bundle(merged_run_root, bundle_spec; run_id = String(run_id), overwrite = true)
    write_posterior_hdf5(paths.posterior_path, samples, parameter_names, first_spec, transport.labels, pathology.timepoints; chain_ids = chain_ids)

    pred_df = DataFrame(region = transport.labels)
    for (j, t) in enumerate(pathology.timepoints)
        pred_df[!, string(t)] = pred_observed[:, j]
    end
    CSV.write(paths.predictions_train_path, pred_df)
    CSV.write(paths.posterior_summary_path, posterior_summary_table(samples, parameter_names))

    diagnostics = Dict(
        "status" => "merged",
        "n_input_runs" => length(run_dirs),
        "n_samples_total" => size(samples, 1),
        "source_run_ids" => [basename(run_dir) for run_dir in run_dirs],
    )
    open(paths.diagnostics_path, "w") do io
        print(io, _json(diagnostics))
    end

    return paths
end
