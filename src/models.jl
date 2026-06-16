const MODEL_PARAMETER_NAMES = Dict(
    "DIFF" => ["rho"],
    "DIFF-R" => ["rho", "alpha", "beta"],
    "DIFF-RF" => ["rho", "alpha", "beta", "gamma"],
    "DIFF-RF-REGIONAL" => ["rho", "alpha", "beta", "gamma"],
    "LOCAL-RF" => ["alpha", "beta", "gamma"],
)

const MODEL_REGIONAL_PARAMETER_NAMES = Dict(
    "DIFF" => String[],
    "DIFF-R" => ["beta"],
    "DIFF-RF" => ["beta", "gamma"],
    "DIFF-RF-REGIONAL" => ["alpha", "beta", "gamma"],
    "LOCAL-RF" => ["beta", "gamma"],
)

const MODEL_STATE_DIMENSIONS = Dict(
    "DIFF" => N -> N,
    "DIFF-R" => N -> N,
    "DIFF-RF" => N -> 2 * N,
    "DIFF-RF-REGIONAL" => N -> 2 * N,
    "LOCAL-RF" => N -> 2 * N,
)

function diff_model!(du, u, p, t; L)
    rho = p[1]
    du .= -rho * L * u
    return nothing
end

function diff_r_model!(du, u, p, t; L)
    N = length(u)
    rho = p[1]
    alpha = p[2]
    beta = @view p[3:(N + 2)]
    du .= -rho * L * u .+ alpha .* u .* (beta .- u)
    return nothing
end

function diff_rf_model!(du, u, p, t; L)
    N = size(L, 1)
    rho = p[1]
    alpha = p[2]
    beta = @view p[3:(N + 2)]
    gamma = @view p[(N + 3):(2 * N + 2)]

    x = @view u[1:N]
    y = @view u[(N + 1):(2 * N)]

    du[1:N] .= -rho * L * x .+ alpha .* x .* (beta .- y .- x)
    du[(N + 1):(2 * N)] .= gamma .* x
    return nothing
end

function diff_rf_regional_model!(du, u, p, t; L)
    N = size(L, 1)
    rho = p[1]
    alpha = @view p[2:(N + 1)]
    beta = @view p[(N + 2):(2 * N + 1)]
    gamma = @view p[(2 * N + 2):(3 * N + 1)]

    x = @view u[1:N]
    y = @view u[(N + 1):(2 * N)]

    du[1:N] .= -rho * L * x .+ alpha .* x .* (beta .- y .- x)
    du[(N + 1):(2 * N)] .= gamma .* x
    return nothing
end

function local_rf_model!(du, u, p, t; L)
    N = div(length(u), 2)
    alpha = p[1]
    beta = @view p[2:(N + 1)]
    gamma = @view p[(N + 2):(2 * N + 1)]

    x = @view u[1:N]
    y = @view u[(N + 1):(2 * N)]

    du[1:N] .= alpha .* x .* (beta .- y .- x)
    du[(N + 1):(2 * N)] .= gamma .* x
    return nothing
end

const MODEL_FUNCTIONS = Dict(
    "DIFF" => diff_model!,
    "DIFF-R" => diff_r_model!,
    "DIFF-RF" => diff_rf_model!,
    "DIFF-RF-REGIONAL" => diff_rf_regional_model!,
    "LOCAL-RF" => local_rf_model!,
)

function initial_state(model_name::AbstractString, N::Integer)
    state_dim = MODEL_STATE_DIMENSIONS[String(model_name)](N)
    return zeros(Float64, state_dim)
end

function default_parameter_vector(model_name::AbstractString, N::Integer)
    name = String(model_name)
    if name == "DIFF"
        return [0.2]
    elseif name == "DIFF-R"
        return vcat([0.1, 0.4], fill(1.0, N))
    elseif name == "DIFF-RF"
        return vcat([0.1, 0.4], fill(1.0, N), fill(0.05, N))
    elseif name == "DIFF-RF-REGIONAL"
        return vcat([0.1], fill(0.4, N), fill(1.0, N), fill(0.05, N))
    elseif name == "LOCAL-RF"
        return vcat([0.4], fill(1.0, N), fill(0.05, N))
    else
        error("Unknown model name: $model_name")
    end
end
