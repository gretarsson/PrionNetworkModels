function make_ode_problem(spec::RunSpec, L::AbstractMatrix, labels::Vector{String}, timepoints::AbstractVector{<:Real}; u0::Union{Nothing,AbstractVector}=nothing)
    N = size(L, 1)
    state0 = isnothing(u0) ? initial_state(spec.model.name, N) : copy(u0)
    model_fn = MODEL_FUNCTIONS[spec.model.name]
    rhs = (du, u, p, t) -> model_fn(du, u, p, t; L=L)
    return ODEProblem(rhs, state0, (0.0, Float64(timepoints[end])))
end

function initial_conditions_for_spec(spec::RunSpec, N::Integer; seed_value=1.0, local_u0=nothing)
    if spec.model.name == "LOCAL-RF"
        values = isnothing(local_u0) ? zeros(Float64, N) : collect(local_u0)
        length(values) == N || error("Expected $N local initial conditions, got $(length(values))")
        T = promote_type(eltype(values), Float64)
        u0 = zeros(T, length(initial_state(spec.model.name, N)))
        u0[1:N] .= values
        return u0
    end

    seed_values = seed_value isa AbstractVector ? collect(seed_value) : fill(seed_value, length(spec.seeding.seed_indices))
    length(seed_values) == length(spec.seeding.seed_indices) || error("Expected $(length(spec.seeding.seed_indices)) seed values, got $(length(seed_values))")

    T = promote_type(eltype(seed_values), Float64)
    u0 = zeros(T, length(initial_state(spec.model.name, N)))
    for (idx, value) in zip(spec.seeding.seed_indices, seed_values)
        u0[idx] = value
    end
    return u0
end

function simulate_trajectory(spec::RunSpec, L::AbstractMatrix, labels::Vector{String}, timepoints::AbstractVector{<:Real}, params::AbstractVector{<:Real}; seed_value=1.0)
    N = size(L, 1)
    seed_values = seed_value isa AbstractVector ? collect(seed_value) : [seed_value]
    T = promote_type(eltype(params), eltype(seed_values), eltype(L), Float64)
    u0 = spec.model.name == "LOCAL-RF" ?
        initial_conditions_for_spec(spec, N; local_u0 = seed_value) :
        initial_conditions_for_spec(spec, N; seed_value = seed_value)
    prob = make_ode_problem(spec, Matrix{T}(L), labels, timepoints; u0=u0)
    sol = solve(prob, Tsit5(); p=collect(params), saveat=Float64.(timepoints), abstol=1e-8, reltol=1e-8)
    return Array(sol)
end
