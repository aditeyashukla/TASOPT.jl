using JuMP
import Ipopt

export IpoptFdOptions, optimize_aircraft_ipopt

"""
    IpoptFdOptions

Controls the IPOPT solve that uses central finite-difference sensitivities with
relative perturbations.
"""
@kwdef struct IpoptFdOptions
    relative_step::Float64 = 1e-3
    min_absolute_step::Float64 = 1e-6
    max_step_reductions::Int = 8
    respect_variable_bounds_in_fd::Bool = false
    tol::Float64 = 1e-6
    dual_inf_tol::Float64 = 1.0
    acceptable_tol::Float64 = 1e-4
    constr_viol_tol::Float64 = 1e-4
    compl_inf_tol::Float64 = 1e-4
    acceptable_iter::Int = 15
    max_iter::Int = 300
    print_level::Int = 5
    progress::Bool = false
    progress_frequency::Int = 10
    progress_print_design_variables::Bool = false
    hessian_approximation::String = "limited-memory"
    derivative_test::String = "none"
    derivative_test_perturbation::Float64 = 1e-6
    derivative_test_tol::Float64 = 1e-4
    derivative_test_print_all::Bool = false
end

mutable struct IpoptFdState
    evaluator::AircraftOptimizationEvaluator
    eval_cache::Dict{Tuple{Vararg{Float64}}, AircraftOptimizationEvaluation}
    fd_x::Union{Nothing, Vector{Float64}}
    fd_objective_gradient::Union{Nothing, Vector{Float64}}
    fd_constraint_jacobian::Union{Nothing, Matrix{Float64}}
end

function IpoptFdState(
    problem::AircraftOptimizationProblem,
    options::IpoptFdOptions,
)
    evaluator = AircraftOptimizationEvaluator(
        problem;
        progress=options.progress,
        progress_frequency=options.progress_frequency,
        progress_print_design_variables=options.progress_print_design_variables,
    )
    return IpoptFdState(
        evaluator,
        Dict{Tuple{Vararg{Float64}}, AircraftOptimizationEvaluation}(),
        nothing,
        nothing,
        nothing,
    )
end

function validate_ipopt_fd_options(options::IpoptFdOptions)
    0.0 < options.relative_step < 1.0 || error("relative_step must lie between 0 and 1.")
    options.min_absolute_step > 0.0 || error("min_absolute_step must be positive.")
    options.max_step_reductions >= 0 || error("max_step_reductions must be nonnegative.")
    options.tol > 0.0 || error("tol must be positive.")
    options.dual_inf_tol > 0.0 || error("dual_inf_tol must be positive.")
    options.acceptable_tol > 0.0 || error("acceptable_tol must be positive.")
    options.constr_viol_tol > 0.0 || error("constr_viol_tol must be positive.")
    options.compl_inf_tol > 0.0 || error("compl_inf_tol must be positive.")
    options.acceptable_iter >= 0 || error("acceptable_iter must be nonnegative.")
    options.max_iter > 0 || error("max_iter must be positive.")
    options.progress_frequency > 0 || error("progress_frequency must be positive.")
    options.derivative_test in ("none", "first-order", "second-order", "only-second-order") ||
        error(
            "derivative_test must be one of " *
            "\"none\", \"first-order\", \"second-order\", or \"only-second-order\".",
        )
    options.derivative_test_perturbation > 0.0 ||
        error("derivative_test_perturbation must be positive.")
    options.derivative_test_tol > 0.0 || error("derivative_test_tol must be positive.")
    return nothing
end

vector_key(x::AbstractVector{<:Real}) = Tuple(Float64.(collect(x)))

function interior_start_value(var::OptimizationDesignVariable)
    span = var.upper - var.lower
    margin = max(1e-6 * span, eps(Float64) * max(abs(var.lower), abs(var.upper), 1.0))
    if margin >= 0.5 * span
        return 0.5 * (var.lower + var.upper)
    end
    return clamp(var.initial, var.lower + margin, var.upper - margin)
end

function interior_start_vector(problem::AircraftOptimizationProblem)
    return [interior_start_value(var) for var in problem.design_variables]
end

function interior_scaled_start_vector(problem::AircraftOptimizationProblem)
    return encode_design_vector(problem, interior_start_vector(problem))
end

function maybe_store_best_model!(
    evaluator::AircraftOptimizationEvaluator,
    evaluation::AircraftOptimizationEvaluation,
    model::Union{Nothing, aircraft},
)
    if is_better_evaluation(evaluation, evaluator.best_candidate)
        evaluator.best_candidate = evaluation
        evaluator.best_candidate_model = model
    end

    if evaluation.feasible && is_better_evaluation(evaluation, evaluator.best_feasible)
        evaluator.best_feasible = evaluation
        evaluator.best_feasible_model = model
    end

    return nothing
end

function evaluate_physical_cached!(
    state::IpoptFdState,
    x::AbstractVector{<:Real},
)
    x_physical = Float64.(collect(x))
    key = vector_key(x_physical)
    if haskey(state.eval_cache, key)
        return state.eval_cache[key]
    end

    evaluation, model = evaluate_aircraft_problem(
        state.evaluator.problem,
        x_physical;
        scaled=false,
        return_model=true,
    )

    push!(state.evaluator.history, evaluation)
    state.evaluator.last_x_scaled = nothing
    state.evaluator.last_evaluation = evaluation
    state.evaluator.last_model = model
    state.eval_cache[key] = evaluation

    maybe_store_best_model!(state.evaluator, evaluation, model)

    progress_row = build_progress_row(state.evaluator, evaluation)
    push!(state.evaluator.progress_rows, progress_row)
    maybe_log_progress!(state.evaluator, evaluation, progress_row)

    return evaluation
end

function evaluate_scaled_cached!(
    state::IpoptFdState,
    x_scaled::AbstractVector{<:Real},
)
    x_scaled_f = Float64.(collect(x_scaled))
    x_physical = decode_design_vector(state.evaluator.problem, x_scaled_f)
    return evaluate_physical_cached!(state, x_physical)
end

function finite_difference_perturbations(
    problem::AircraftOptimizationProblem,
    x::Vector{Float64},
    index::Int,
    relative_step::Float64,
    options::IpoptFdOptions,
)
    var = problem.design_variables[index]
    x_center = x[index]
    if abs(x_center) >= options.min_absolute_step
        x_plus_value = x_center * (1.0 + relative_step)
        x_minus_value = x_center * (1.0 - relative_step)
    else
        x_plus_value = x_center + options.min_absolute_step
        x_minus_value = x_center - options.min_absolute_step
    end

    if options.respect_variable_bounds_in_fd
        x_plus_value = clamp(x_plus_value, var.lower, var.upper)
        x_minus_value = clamp(x_minus_value, var.lower, var.upper)
    end

    if x_plus_value == x_minus_value
        x_plus_value = x_center + options.min_absolute_step
        x_minus_value = x_center - options.min_absolute_step
        if options.respect_variable_bounds_in_fd
            x_plus_value = clamp(x_plus_value, var.lower, var.upper)
            x_minus_value = clamp(x_minus_value, var.lower, var.upper)
        end
    end

    x_plus_value != x_minus_value || error(
        "Central finite-difference step for $(var.name) collapsed at x=$(x_center).",
    )

    x_plus = copy(x)
    x_minus = copy(x)
    x_plus[index] = x_plus_value
    x_minus[index] = x_minus_value
    return x_plus, x_minus
end

function evaluate_finite_difference_pair!(
    state::IpoptFdState,
    x_physical::Vector{Float64},
    index::Int,
    options::IpoptFdOptions,
)
    problem = state.evaluator.problem
    relative_step = options.relative_step
    last_trial = nothing

    for _ in 0:options.max_step_reductions
        x_plus, x_minus = finite_difference_perturbations(
            problem,
            x_physical,
            index,
            relative_step,
            options,
        )
        plus_evaluation = evaluate_physical_cached!(state, x_plus)
        minus_evaluation = evaluate_physical_cached!(state, x_minus)
        last_trial = (x_plus, x_minus, plus_evaluation, minus_evaluation)

        if plus_evaluation.status == :ok && minus_evaluation.status == :ok
            return last_trial
        end

        relative_step *= 0.5
    end

    return last_trial
end

function ensure_fd_cache!(
    state::IpoptFdState,
    x_scaled::AbstractVector{<:Real},
    options::IpoptFdOptions,
)
    x_scaled_f = Float64.(collect(x_scaled))
    if !isnothing(state.fd_x) && state.fd_x == x_scaled_f
        return nothing
    end

    problem = state.evaluator.problem
    x_physical = decode_design_vector(problem, x_scaled_f)
    center_evaluation = evaluate_physical_cached!(state, x_physical)
    nvars = length(problem.design_variables)
    nconstraints = length(problem.constraints)

    objective_gradient = zeros(nvars)
    constraint_jacobian = zeros(nconstraints, nvars)

    for j in 1:nvars
        x_plus, x_minus, plus_evaluation, minus_evaluation = evaluate_finite_difference_pair!(
            state,
            x_physical,
            j,
            options,
        )

        denominator = x_plus[j] - x_minus[j]
        physical_gradient =
            (plus_evaluation.objective - minus_evaluation.objective) / denominator
        variable_scale = problem.design_variables[j].upper - problem.design_variables[j].lower
        objective_gradient[j] = physical_gradient * variable_scale

        if nconstraints > 0
            constraint_jacobian[:, j] .=
                ((plus_evaluation.constraint_values .- minus_evaluation.constraint_values) ./ denominator) .* variable_scale
        end
    end

    state.fd_x = x_scaled_f
    state.fd_objective_gradient = objective_gradient
    state.fd_constraint_jacobian = constraint_jacobian

    return center_evaluation
end

function objective_operator_function(state::IpoptFdState)
    return function (x...)
        evaluation = evaluate_scaled_cached!(state, Float64.(collect(x)))
        return evaluation.objective
    end
end

function objective_gradient_function(
    state::IpoptFdState,
    options::IpoptFdOptions,
)
    return function (gradient::AbstractVector, x...)
        ensure_fd_cache!(state, Float64.(collect(x)), options)
        copyto!(gradient, state.fd_objective_gradient)
        return nothing
    end
end

function constraint_operator_function(
    state::IpoptFdState,
    constraint_index::Int,
)
    return function (x...)
        evaluation = evaluate_scaled_cached!(state, Float64.(collect(x)))
        return evaluation.constraint_values[constraint_index]
    end
end

function constraint_gradient_function(
    state::IpoptFdState,
    options::IpoptFdOptions,
    constraint_index::Int,
)
    return function (gradient::AbstractVector, x...)
        ensure_fd_cache!(state, Float64.(collect(x)), options)
        copyto!(gradient, @view(state.fd_constraint_jacobian[constraint_index, :]))
        return nothing
    end
end

function ipopt_attributes!(
    model::JuMP.Model,
    options::IpoptFdOptions,
)
    JuMP.set_attribute(model, "tol", options.tol)
    JuMP.set_attribute(model, "dual_inf_tol", options.dual_inf_tol)
    JuMP.set_attribute(model, "acceptable_tol", options.acceptable_tol)
    JuMP.set_attribute(model, "constr_viol_tol", options.constr_viol_tol)
    JuMP.set_attribute(model, "compl_inf_tol", options.compl_inf_tol)
    JuMP.set_attribute(model, "acceptable_iter", options.acceptable_iter)
    JuMP.set_attribute(model, "max_iter", options.max_iter)
    JuMP.set_attribute(model, "print_level", options.print_level)
    JuMP.set_attribute(model, "hessian_approximation", options.hessian_approximation)
    if options.derivative_test != "none"
        JuMP.set_attribute(model, "derivative_test", options.derivative_test)
        JuMP.set_attribute(
            model,
            "derivative_test_perturbation",
            options.derivative_test_perturbation,
        )
        JuMP.set_attribute(model, "derivative_test_tol", options.derivative_test_tol)
        JuMP.set_attribute(
            model,
            "derivative_test_print_all",
            options.derivative_test_print_all ? "yes" : "no",
        )
    end
    if options.print_level <= 0
        JuMP.set_attribute(model, "sb", "yes")
    end
    return nothing
end

"""
    optimize_aircraft_ipopt(problem; options=IpoptFdOptions())

Runs an IPOPT solve in normalized `[0, 1]` design-variable space using central
finite-difference sensitivities with relative perturbations.
"""
function optimize_aircraft_ipopt(
    problem::AircraftOptimizationProblem;
    options::IpoptFdOptions = IpoptFdOptions(),
)
    validate_ipopt_fd_options(options)

    state = IpoptFdState(problem, options)
    start_values = interior_scaled_start_vector(problem)

    # Prime the cache so the result always contains at least one evaluation.
    evaluate_scaled_cached!(state, start_values)

    model = JuMP.Model(Ipopt.Optimizer)
    ipopt_attributes!(model, options)

    nvars = length(problem.design_variables)

    @variable(
        model,
        x[i = 1:nvars],
        lower_bound = 0.0,
        upper_bound = 1.0,
        start = start_values[i],
    )

    objective_operator = JuMP.add_nonlinear_operator(
        model,
        nvars,
        objective_operator_function(state),
        objective_gradient_function(state, options);
        name=:aircraft_fd_objective,
    )
    @objective(model, Min, objective_operator(x...))

    for constraint_index in eachindex(problem.constraints)
        operator_name = Symbol("aircraft_fd_constraint_", constraint_index)
        constraint_operator = JuMP.add_nonlinear_operator(
            model,
            nvars,
            constraint_operator_function(state, constraint_index),
            constraint_gradient_function(state, options, constraint_index);
            name=operator_name,
        )
        @constraint(model, constraint_operator(x...) <= 0.0)
    end

    JuMP.optimize!(model)

    if JuMP.has_values(model)
        evaluate_scaled_cached!(state, JuMP.value.(x))
    end

    best_candidate = state.evaluator.best_candidate
    best_candidate_model = state.evaluator.best_candidate_model
    best_feasible = state.evaluator.best_feasible
    best_feasible_model = state.evaluator.best_feasible_model

    optimum = isnothing(best_feasible) ? best_candidate : best_feasible
    optimum_model = isnothing(best_feasible) ? best_candidate_model : best_feasible_model

    return AircraftOptimizationResult(
        :IPOPT_FD,
        Any[(
            termination_status=JuMP.termination_status(model),
            primal_status=JuMP.primal_status(model),
            raw_status=JuMP.raw_status(model),
        )],
        optimum,
        optimum_model,
        best_feasible,
        best_feasible_model,
        state.evaluator.history,
        state.evaluator.progress_rows,
        length(state.evaluator.history),
    )
end
