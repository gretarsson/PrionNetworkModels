function make_ode_problem(spec::RunSpec, L::AbstractMatrix, labels::Vector{String}, timepoints::AbstractVector{<:Real}; u0::Union{Nothing,AbstractVector}=nothing)
    N = size(L, 1)
    state0 = isnothing(u0) ? initial_state(spec.model.name, N) : copy(u0)
    model_fn = MODEL_FUNCTIONS[spec.model.name]
    rhs = (du, u, p, t) -> model_fn(du, u, p, t; L=L)
    return ODEProblem(rhs, state0, (0.0, Float64(timepoints[end])))
end

function initial_conditions_for_spec(spec::RunSpec, N::Integer; seed_value::Real=1.0)
    T = promote_type(typeof(seed_value), Float64)
    u0 = zeros(T, length(initial_state(spec.model.name, N)))
    for idx in spec.seeding.seed_indices
        u0[idx] = seed_value
    end
    return u0
end

function simulate_trajectory(spec::RunSpec, L::AbstractMatrix, labels::Vector{String}, timepoints::AbstractVector{<:Real}, params::AbstractVector{<:Real}; seed_value::Real=1.0)
    N = size(L, 1)
    T = promote_type(eltype(params), typeof(seed_value), eltype(L), Float64)
    u0 = initial_conditions_for_spec(spec, N; seed_value=seed_value)
    prob = make_ode_problem(spec, Matrix{T}(L), labels, timepoints; u0=u0)
    sol = solve(prob, Tsit5(); p=collect(params), saveat=Float64.(timepoints), abstol=1e-8, reltol=1e-8)
    return Array(sol)
end
