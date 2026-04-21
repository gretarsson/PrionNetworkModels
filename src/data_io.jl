function process_pathology(path_csv::AbstractString; network_csv::Union{Nothing,String}=nothing)
    df = CSV.read(path_csv, DataFrame; missingstring=["NA", ""])

    time_col = names(df)[2]
    region_order = names(df)[3:end]

    if network_csv !== nothing
        ndf = CSV.read(network_csv, DataFrame)
        network_regions = names(ndf)[2:end]
        missing_regions = setdiff(network_regions, region_order)
        if !isempty(missing_regions)
            for region in missing_regions
                df[!, Symbol(region)] = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
            end
        end
        region_order = network_regions
    end

    timepoints = sort(unique(df[!, time_col]))
    groups = groupby(df, time_col)
    max_samples = maximum(nrow(g) for g in groups)

    n_regions = length(region_order)
    n_times = length(timepoints)
    data = Array{Union{Float64, Missing}}(undef, n_regions, n_times, max_samples)
    fill!(data, missing)

    group_map = Dict(g[1, time_col] => g for g in groups)
    for (j, timepoint) in enumerate(timepoints)
        sub = group_map[timepoint]
        for s in 1:nrow(sub)
            for (k, region) in enumerate(region_order)
                data[k, j, s] = sub[s, region]
            end
        end
    end

    return (
        data = data,
        timepoints = Float64.(timepoints),
        labels = String.(region_order),
        max_samples = max_samples,
    )
end

function mean_over_replicates(data::AbstractArray)
    if ndims(data) == 2
        return data
    end
    R, T, _ = size(data)
    out = Array{Union{Float64, Missing}}(missing, R, T)
    for r in 1:R, t in 1:T
        out[r, t] = mean(skipmissing(data[r, t, :]))
    end
    out[isnan.(out)] .= missing
    return out
end

function summarize_over_replicates(data::AbstractArray)
    mean_mat = mean_over_replicates(data)

    if ndims(data) == 2
        counts = Float64.(isfinite.(Float64.(coalesce.(data, NaN))))
        zeros_mat = zeros(Float64, size(mean_mat))
        return (
            mean = Float64.(coalesce.(mean_mat, NaN)),
            sd = zeros_mat,
            se = zeros_mat,
            n = counts,
        )
    end

    R, T, _ = size(data)
    sd = fill(NaN, R, T)
    se = fill(NaN, R, T)
    n = zeros(Float64, R, T)

    for r in 1:R, t in 1:T
        values = collect(skipmissing(data[r, t, :]))
        n[r, t] = length(values)
        if isempty(values)
            continue
        elseif length(values) == 1
            sd[r, t] = 0.0
            se[r, t] = 0.0
        else
            sd[r, t] = std(values)
            se[r, t] = sd[r, t] / sqrt(length(values))
        end
    end

    return (
        mean = Float64.(coalesce.(mean_mat, NaN)),
        sd = sd,
        se = se,
        n = n,
    )
end
