function parameter_count(model_name::AbstractString, N::Integer; parameter_sharing::AbstractString="independent")
    sharing = String(parameter_sharing)
    groups = sharing == "bilateral_pairs" ? div(N, 2) : N
    name = String(model_name)
    if name == "DIFF"
        return 1
    elseif name == "DIFF-R"
        return 2 + groups
    elseif name == "DIFF-RF"
        return 2 + 2 * groups
    elseif name == "DIFF-RF-REGIONAL"
        return 1 + 3 * groups
    elseif name == "LOCAL-RF"
        return 1 + 3 * groups
    else
        error("Unknown model name: $model_name")
    end
end
