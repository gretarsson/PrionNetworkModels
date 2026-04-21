function observed_matrix(data::AbstractArray)
    out = mean_over_replicates(data)
    return Float64.(coalesce.(out, NaN))
end
