using NLopt
using Random

export OptimizationDesignVariable, OptimizationConstraint,
       MissionOptimizationSettings, AircraftOptimizationProblem,
       AircraftOptimizationEvaluation, AircraftOptimizationResult,
       CobylaOptions, build_weighted_pfei_problem,
       initial_physical_design_vector, initial_scaled_design_vector,
       encode_design_vector, decode_design_vector,
       evaluate_aircraft_problem, optimize_aircraft,
       plot_optimization_history,
       weighted_pfei_objective, balanced_field_length_limit,
       fuel_capacity_constraints, balanced_field_length_constraints,
       wing_span_constraint, top_of_climb_gradient_constraint,
       tt3_constraint, vane_temperature_constraint,
       takeoff_weight_constraint, fan_diameter_constraint,
       fan_pr_less_than_lpc_pr_constraints,
       standard_weighted_pfei_constraints,
       aspect_ratio_variable, wing_sweep_variable,
       inboard_taper_ratio_variable, outboard_taper_ratio_variable,
       root_thickness_to_chord_variable, spanbreak_thickness_to_chord_variable,
       mission_cruise_altitude_variable, mission_cruise_cl_variable,
       mission_rcls_variable, mission_rclt_variable,
       shared_mission_rcls_variable, shared_mission_rclt_variable,
       shared_cruise_turbine_temperature_variable,
       shared_takeoff_turbine_temperature_variable,
       shared_fan_pr_variable, shared_lpc_pr_variable,
       shared_hpc_pr_variable, shared_bpr_variable,
       deviation_height_variable, deviation_length_variable,
       deviation_start_variable

"""
    OptimizationDesignVariable(name, lower, upper, initial, setter!)

Defines one scalar design variable for an aircraft optimization problem.

`setter!` must have signature `(ac::aircraft, value::Float64) -> nothing`.
"""
struct OptimizationDesignVariable
    name::String
    lower::Float64
    upper::Float64
    initial::Float64
    setter!::Function
end

function OptimizationDesignVariable(
    name::AbstractString,
    lower::Real,
    upper::Real,
    initial::Real,
    setter!::Function,
)
    lower_f = Float64(lower)
    upper_f = Float64(upper)
    initial_f = Float64(initial)
    isfinite(lower_f) || error("Lower bound for $name must be finite.")
    isfinite(upper_f) || error("Upper bound for $name must be finite.")
    lower_f < upper_f || error("Lower bound must be smaller than upper bound for $name.")
    lower_f <= initial_f <= upper_f || error("Initial value must lie inside bounds for $name.")
    return OptimizationDesignVariable(String(name), lower_f, upper_f, initial_f, setter!)
end

"""
    OptimizationConstraint(name, func; tolerance=1e-3)

Defines one nonlinear inequality constraint `g(x) <= 0`.

`func` must have signature `(ac::aircraft) -> Float64`.
"""
struct OptimizationConstraint
    name::String
    func::Function
    tolerance::Float64
end

function OptimizationConstraint(
    name::AbstractString,
    func::Function;
    tolerance::Real = 1e-3,
)
    tolerance_f = Float64(tolerance)
    tolerance_f >= 0.0 || error("Constraint tolerance for $name must be nonnegative.")
    return OptimizationConstraint(String(name), func, tolerance_f)
end

"""
    MissionOptimizationSettings

Controls how the design mission is sized and how off-design missions are run.
"""
@kwdef struct MissionOptimizationSettings
    design_iter::Int = 50
    offdesign_iter::Int = 50
    printiter::Bool = false
    printTO::Bool = false
    initializes_engine::Bool = true
    mission_modes::Dict{Int, Symbol} = Dict{Int, Symbol}()
    offdesign_missions::Vector{Int} = Int[]
    parallel_offdesign::Bool = false
    parallel_offdesign_workers::Int = 0
end

"""
    AircraftOptimizationProblem

Stores the base aircraft, decision variables, objective, and constraints.
"""
struct AircraftOptimizationProblem
    base_model::aircraft
    design_variables::Vector{OptimizationDesignVariable}
    objective::Function
    constraints::Vector{OptimizationConstraint}
    run_settings::MissionOptimizationSettings
    failure_objective::Float64
    failure_constraint::Float64
    progress_row_builder::Union{Nothing, Function}
end

function AircraftOptimizationProblem(
    base_model::aircraft,
    design_variables::Vector{OptimizationDesignVariable},
    objective::Function;
    constraints::Vector{OptimizationConstraint} = OptimizationConstraint[],
    run_settings::MissionOptimizationSettings = MissionOptimizationSettings(),
    failure_objective::Real = 1e12,
    failure_constraint::Real = 1e6,
    progress_row_builder::Union{Nothing, Function} = nothing,
)
    isempty(design_variables) && error("AircraftOptimizationProblem requires at least one design variable.")
    validate_mission_settings(base_model, run_settings)
    return AircraftOptimizationProblem(
        base_model,
        design_variables,
        objective,
        constraints,
        run_settings,
        Float64(failure_objective),
        Float64(failure_constraint),
        progress_row_builder,
    )
end

"""
    AircraftOptimizationEvaluation

One cached aircraft evaluation at a point in the scaled design space.
"""
struct AircraftOptimizationEvaluation
    x_scaled::Vector{Float64}
    x_physical::Vector{Float64}
    objective::Float64
    constraint_values::Vector{Float64}
    feasible::Bool
    max_violation::Float64
    status::Symbol
    message::String
end

"""
    CobylaOptions

Controls the NLopt derivative-free solve.
"""
@kwdef struct CobylaOptions
    algorithm::Symbol = :LN_COBYLA
    local_algorithm::Union{Nothing, Symbol} = nothing
    initial_step::Float64 = 0.05
    xtol_rel::Float64 = 1e-3
    ftol_rel::Float64 = 1e-6
    maxeval::Int = 1500
    num_starts::Int = 1
    start_points_physical::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    start_points_scaled::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    seed::Int = 1234
    progress::Bool = false
    progress_frequency::Int = 10
    progress_print_design_variables::Bool = false
end

"""
    AircraftOptimizationResult

Stores the best candidate found by `optimize_aircraft`.
"""
struct AircraftOptimizationResult
    algorithm::Symbol
    start_return_codes::Vector{Any}
    optimum::AircraftOptimizationEvaluation
    optimum_model::Union{Nothing, aircraft}
    best_feasible::Union{Nothing, AircraftOptimizationEvaluation}
    best_feasible_model::Union{Nothing, aircraft}
    history::Vector{AircraftOptimizationEvaluation}
    progress_rows::Vector{Vector{Pair{String, Any}}}
    num_evaluations::Int
end

mutable struct AircraftOptimizationEvaluator
    problem::AircraftOptimizationProblem
    progress::Bool
    progress_frequency::Int
    progress_print_design_variables::Bool
    history::Vector{AircraftOptimizationEvaluation}
    last_x_scaled::Union{Nothing, Vector{Float64}}
    last_evaluation::Union{Nothing, AircraftOptimizationEvaluation}
    last_model::Union{Nothing, aircraft}
    best_candidate::Union{Nothing, AircraftOptimizationEvaluation}
    best_candidate_model::Union{Nothing, aircraft}
    best_feasible::Union{Nothing, AircraftOptimizationEvaluation}
    best_feasible_model::Union{Nothing, aircraft}
    progress_rows::Vector{Vector{Pair{String, Any}}}
end

function AircraftOptimizationEvaluator(
    problem::AircraftOptimizationProblem;
    progress::Bool = false,
    progress_frequency::Int = 10,
    progress_print_design_variables::Bool = false,
)
    return AircraftOptimizationEvaluator(
        problem,
        progress,
        progress_frequency,
        progress_print_design_variables,
        AircraftOptimizationEvaluation[],
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        Vector{Pair{String, Any}}[],
    )
end

function validate_mission_settings(base_model::aircraft, settings::MissionOptimizationSettings)
    nmissions = size(base_model.parm, 2)
    settings.design_iter > 1 || error("design_iter must be greater than 1.")
    settings.offdesign_iter > 1 || error("offdesign_iter must be greater than 1.")

    for (mission, mode) in settings.mission_modes
        2 <= mission <= nmissions || error("Mission mode specified for invalid off-design mission $mission.")
        mode in (:altitude, :CL) || error("Mission mode for mission $mission must be :altitude or :CL.")
    end

    for mission in settings.offdesign_missions
        2 <= mission <= nmissions || error("Off-design mission index $mission is out of bounds.")
    end

    return nothing
end

function encode_value(var::OptimizationDesignVariable, x::Real)
    return (Float64(x) - var.lower) / (var.upper - var.lower)
end

function decode_value(var::OptimizationDesignVariable, x_scaled::Real)
    return var.lower + Float64(x_scaled) * (var.upper - var.lower)
end

"""
    initial_physical_design_vector(problem)

Returns the physical initial values for each design variable.
"""
function initial_physical_design_vector(problem::AircraftOptimizationProblem)
    return [var.initial for var in problem.design_variables]
end

"""
    initial_scaled_design_vector(problem)

Returns the initial point scaled to `[0, 1]`.
"""
function initial_scaled_design_vector(problem::AircraftOptimizationProblem)
    return encode_design_vector(problem, initial_physical_design_vector(problem))
end

"""
    encode_design_vector(problem, x_physical)

Maps a physical design vector into `[0, 1]^n`.
"""
function encode_design_vector(problem::AircraftOptimizationProblem, x_physical::AbstractVector{<:Real})
    length(x_physical) == length(problem.design_variables) || error("Physical design vector has the wrong length.")
    return [encode_value(problem.design_variables[i], x_physical[i]) for i in eachindex(problem.design_variables)]
end

"""
    decode_design_vector(problem, x_scaled)

Maps a scaled design vector from `[0, 1]^n` into physical units.
"""
function decode_design_vector(problem::AircraftOptimizationProblem, x_scaled::AbstractVector{<:Real})
    length(x_scaled) == length(problem.design_variables) || error("Scaled design vector has the wrong length.")
    return [decode_value(problem.design_variables[i], x_scaled[i]) for i in eachindex(problem.design_variables)]
end

function offdesign_missions(problem::AircraftOptimizationProblem)
    if isempty(problem.run_settings.offdesign_missions)
        nmissions = size(problem.base_model.parm, 2)
        return collect(2:nmissions)
    end
    return problem.run_settings.offdesign_missions
end

function mission_mode_string(problem::AircraftOptimizationProblem, mission::Int)
    mode = get(problem.run_settings.mission_modes, mission, :altitude)
    return mode == :CL ? "CL" : "altitude"
end

function apply_design_variables!(
    ac::aircraft,
    design_variables::Vector{OptimizationDesignVariable},
    x_physical::Vector{Float64},
)
    for (var, value) in zip(design_variables, x_physical)
        var.setter!(ac, value)
    end
    return nothing
end

function build_success_evaluation(
    problem::AircraftOptimizationProblem,
    x_scaled::Vector{Float64},
    x_physical::Vector{Float64},
    objective::Float64,
    constraint_values::Vector{Float64},
)
    feasibility_thresholds = [constraint.tolerance for constraint in problem.constraints]
    feasible = all(
        constraint_values[i] <= feasibility_thresholds[i]
        for i in eachindex(constraint_values)
    )
    max_violation = isempty(constraint_values) ? 0.0 : maximum(max.(0.0, constraint_values))
    return AircraftOptimizationEvaluation(
        copy(x_scaled),
        copy(x_physical),
        objective,
        constraint_values,
        feasible,
        max_violation,
        :ok,
        "",
    )
end

function build_failure_evaluation(
    problem::AircraftOptimizationProblem,
    x_scaled::Vector{Float64},
    x_physical::Vector{Float64},
    message::String,
)
    constraint_values = fill(problem.failure_constraint, length(problem.constraints))
    max_violation = isempty(constraint_values) ? 0.0 : maximum(max.(0.0, constraint_values))
    return AircraftOptimizationEvaluation(
        copy(x_scaled),
        copy(x_physical),
        problem.failure_objective,
        constraint_values,
        false,
        max_violation,
        :failed,
        message,
    )
end

function is_better_evaluation(
    candidate::AircraftOptimizationEvaluation,
    incumbent::Union{Nothing, AircraftOptimizationEvaluation},
)
    isnothing(incumbent) && return true

    if candidate.feasible != incumbent.feasible
        return candidate.feasible
    end

    if candidate.feasible
        return candidate.objective < incumbent.objective
    end

    if candidate.max_violation != incumbent.max_violation
        return candidate.max_violation < incumbent.max_violation
    end

    return candidate.objective < incumbent.objective
end

function print_design_variable_progress(
    problem::AircraftOptimizationProblem,
    evaluation::AircraftOptimizationEvaluation,
)
    for (var, value) in zip(problem.design_variables, evaluation.x_physical)
        println(
            "  ",
            rpad(design_variable_display_title(var.name), 36),
            design_variable_display_value(var.name, value),
        )
    end
    return nothing
end

function design_variable_display_unit(name::AbstractString)
    if startswith(name, "cruise_altitude_") || startswith(name, "deviation_height_")
        return "ft"
    elseif startswith(name, "deviation_length_") || startswith(name, "deviation_start_")
        return "nmi"
    end
    return ""
end

function design_variable_display_value(name::AbstractString, value::Real)
    value_f = Float64(value)
    if startswith(name, "cruise_altitude_") || startswith(name, "deviation_height_")
        return value_f / ft_to_m
    elseif startswith(name, "deviation_length_") || startswith(name, "deviation_start_")
        return value_f / nmi_to_m
    end
    return value_f
end

function design_variable_display_name(name::AbstractString)
    unit = design_variable_display_unit(name)
    return isempty(unit) ? String(name) : "$(name)_$(unit)"
end

function design_variable_display_title(name::AbstractString)
    unit = design_variable_display_unit(name)
    return isempty(unit) ? String(name) : "$(name) ($(unit))"
end

function design_variable_display_ylabel(name::AbstractString)
    unit = design_variable_display_unit(name)
    return isempty(unit) ? "Value" : unit
end

function format_progress_value(value)
    if isa(value, Symbol)
        return String(value)
    elseif isa(value, AbstractString)
        return replace(value, '\n' => ' ', ',' => ';')
    elseif isnothing(value)
        return ""
    else
        return string(value)
    end
end

function build_progress_row(
    evaluator::AircraftOptimizationEvaluator,
    evaluation::AircraftOptimizationEvaluation,
)
    row = Pair{String, Any}[
        "eval" => length(evaluator.history),
        "status" => evaluation.status,
        "feasible" => evaluation.feasible,
        "objective" => evaluation.objective,
        "max_violation" => evaluation.max_violation,
    ]

    if !isnothing(evaluator.problem.progress_row_builder)
        extra_row = evaluator.problem.progress_row_builder(
            evaluation,
            evaluator.last_model,
            evaluator.problem,
        )
        append!(row, Pair{String, Any}[String(name) => value for (name, value) in extra_row])
    end

    return row
end

function print_progress_row!(
    evaluator::AircraftOptimizationEvaluator,
    row::Vector{Pair{String, Any}},
)
    if length(evaluator.history) == 1
        println(join(first.(row), ","))
    end
    println(join([format_progress_value(last(pair)) for pair in row], ","))
    return nothing
end

function maybe_log_progress!(
    evaluator::AircraftOptimizationEvaluator,
    evaluation::AircraftOptimizationEvaluation,
    row::Vector{Pair{String, Any}},
)
    if !evaluator.progress
        return nothing
    end

    nevals = length(evaluator.history)
    if nevals == 1 || nevals % evaluator.progress_frequency == 0
        if !isnothing(evaluator.problem.progress_row_builder)
            print_progress_row!(evaluator, row)
        else
            println(
                "Eval $(nevals): f=$(round(evaluation.objective, digits=6)) " *
                "max_g=$(round(evaluation.max_violation, digits=6)) " *
                "status=$(evaluation.status)",
            )
        end

        if isnothing(evaluator.problem.progress_row_builder) && evaluator.progress_print_design_variables
            print_design_variable_progress(evaluator.problem, evaluation)
            println()
        end
    end

    return nothing
end

"""
    evaluate_aircraft_problem(problem, x; scaled=true, return_model=false)

Evaluates the full aircraft optimization problem at `x`.
"""
function evaluate_aircraft_problem(
    problem::AircraftOptimizationProblem,
    x::AbstractVector{<:Real};
    scaled::Bool = true,
    return_model::Bool = false,
)
    x_scaled = scaled ? Float64.(collect(x)) : encode_design_vector(problem, x)
    x_physical = scaled ? decode_design_vector(problem, x_scaled) : Float64.(collect(x))
    ac = deepcopy(problem.base_model)

    try
        apply_design_variables!(ac, problem.design_variables, x_physical)

        size_aircraft!(
            ac;
            iter=problem.run_settings.design_iter,
            printiter=problem.run_settings.printiter,
        )

        missions = offdesign_missions(problem)
        if problem.run_settings.parallel_offdesign && length(missions) > 1
            mode_strings = [mission_mode_string(problem, mission) for mission in missions]
            if all(==(first(mode_strings)), mode_strings)
                fly_mission_parallel!(
                    ac,
                    missions;
                    itermax=problem.run_settings.offdesign_iter,
                    initializes_engine=problem.run_settings.initializes_engine,
                    opt_prescribed_cruise_parameter=first(mode_strings),
                    printTO=problem.run_settings.printTO,
                    nworkers=problem.run_settings.parallel_offdesign_workers > 0 ?
                        problem.run_settings.parallel_offdesign_workers : length(missions),
                    add_workers=true,
                    cleanup_workers=false,
                )
            else
                for mission in missions
                    fly_mission!(
                        ac,
                        mission;
                        itermax=problem.run_settings.offdesign_iter,
                        initializes_engine=problem.run_settings.initializes_engine,
                        opt_prescribed_cruise_parameter=mission_mode_string(problem, mission),
                        printTO=problem.run_settings.printTO,
                    )
                end
            end
        else
            for mission in missions
                fly_mission!(
                    ac,
                    mission;
                    itermax=problem.run_settings.offdesign_iter,
                    initializes_engine=problem.run_settings.initializes_engine,
                    opt_prescribed_cruise_parameter=mission_mode_string(problem, mission),
                    printTO=problem.run_settings.printTO,
                )
            end
        end

        objective = Float64(problem.objective(ac))
        constraint_values = Float64[Float64(constraint.func(ac)) for constraint in problem.constraints]

        if !isfinite(objective) || any(!isfinite, constraint_values)
            message = "Objective or constraint evaluation returned a non-finite value."
            evaluation = build_failure_evaluation(problem, x_scaled, x_physical, message)
            return return_model ? (evaluation, nothing) : evaluation
        end

        evaluation = build_success_evaluation(problem, x_scaled, x_physical, objective, constraint_values)
        return return_model ? (evaluation, ac) : evaluation
    catch err
        message = sprint(showerror, err)
        evaluation = build_failure_evaluation(problem, x_scaled, x_physical, message)
        return return_model ? (evaluation, nothing) : evaluation
    end
end

function evaluate_cached!(
    evaluator::AircraftOptimizationEvaluator,
    x_scaled::AbstractVector{<:Real},
)
    if !isnothing(evaluator.last_x_scaled) && evaluator.last_x_scaled == x_scaled
        return evaluator.last_evaluation
    end

    evaluation, model = evaluate_aircraft_problem(
        evaluator.problem,
        x_scaled;
        scaled=true,
        return_model=true,
    )

    push!(evaluator.history, evaluation)
    evaluator.last_x_scaled = Float64.(collect(x_scaled))
    evaluator.last_evaluation = evaluation
    evaluator.last_model = model

    if is_better_evaluation(evaluation, evaluator.best_candidate)
        evaluator.best_candidate = evaluation
        evaluator.best_candidate_model = model
    end

    if evaluation.feasible && is_better_evaluation(evaluation, evaluator.best_feasible)
        evaluator.best_feasible = evaluation
        evaluator.best_feasible_model = model
    end

    progress_row = build_progress_row(evaluator, evaluation)
    push!(evaluator.progress_rows, progress_row)

    maybe_log_progress!(evaluator, evaluation, progress_row)
    return evaluation
end

function objective_callback!(
    evaluator::AircraftOptimizationEvaluator,
    x::AbstractVector{<:Real},
)
    evaluation = evaluate_cached!(evaluator, x)
    return evaluation.objective
end

function constraint_callback!(
    evaluator::AircraftOptimizationEvaluator,
    x::AbstractVector{<:Real},
    constraint_index::Int,
)
    evaluation = evaluate_cached!(evaluator, x)
    return evaluation.constraint_values[constraint_index]
end

function configure_nlopt_optimizer!(
    opt::NLopt.Opt,
    nvars::Int,
    options::CobylaOptions,
)
    opt.lower_bounds = zeros(nvars)
    opt.upper_bounds = ones(nvars)
    opt.initial_step = fill(options.initial_step, nvars)
    opt.xtol_rel = options.xtol_rel
    opt.ftol_rel = options.ftol_rel
    opt.maxeval = options.maxeval
    return opt
end

function run_nlopt_start(
    problem::AircraftOptimizationProblem,
    start::Vector{Float64},
    options::CobylaOptions,
)
    evaluator = AircraftOptimizationEvaluator(
        problem;
        progress=options.progress,
        progress_frequency=options.progress_frequency,
        progress_print_design_variables=options.progress_print_design_variables,
    )

    nvars = length(problem.design_variables)
    opt = configure_nlopt_optimizer!(
        NLopt.Opt(options.algorithm, nvars),
        nvars,
        options,
    )
    opt.min_objective = (x, grad) -> objective_callback!(evaluator, x)

    if !isnothing(options.local_algorithm)
        local_opt = configure_nlopt_optimizer!(
            NLopt.Opt(options.local_algorithm, nvars),
            nvars,
            options,
        )
        NLopt.local_optimizer!(opt, local_opt)
    end

    for (idx, constraint) in enumerate(problem.constraints)
        NLopt.inequality_constraint!(
            opt,
            (x, grad) -> constraint_callback!(evaluator, x, idx),
            constraint.tolerance,
        )
    end

    _, optx, ret = NLopt.optimize(opt, start)
    final_evaluation = evaluate_cached!(evaluator, optx)

    return (
        return_code = ret,
        evaluation = final_evaluation,
        model = final_evaluation.feasible ? evaluator.best_feasible_model : evaluator.best_candidate_model,
        best_candidate = evaluator.best_candidate,
        best_candidate_model = evaluator.best_candidate_model,
        best_feasible = evaluator.best_feasible,
        best_feasible_model = evaluator.best_feasible_model,
        history = evaluator.history,
        progress_rows = evaluator.progress_rows,
    )
end

function multistart_points(problem::AircraftOptimizationProblem, options::CobylaOptions)
    nvars = length(problem.design_variables)
    starts = Vector{Vector{Float64}}()
    seen = Set{Tuple{Vararg{Int}}}()

    function add_start!(x_scaled)
        length(x_scaled) == nvars || error("COBYLA start point has the wrong length.")
        x = clamp.(Float64.(collect(x_scaled)), 0.0, 1.0)
        key = Tuple(round.(Int, x .* 1.0e10))
        if !(key in seen)
            push!(starts, x)
            push!(seen, key)
        end
        return nothing
    end

    add_start!(initial_scaled_design_vector(problem))

    for x_physical in options.start_points_physical
        add_start!(encode_design_vector(problem, x_physical))
    end

    for x_scaled in options.start_points_scaled
        add_start!(x_scaled)
    end

    if options.num_starts > 1
        rng = MersenneTwister(options.seed)
        for _ = 2:options.num_starts
            add_start!(rand(rng, nvars))
        end
    end

    return starts
end

"""
    optimize_aircraft(problem; options=CobylaOptions())

Runs an NLopt solve on the normalized design space.
"""
function optimize_aircraft(
    problem::AircraftOptimizationProblem;
    options::CobylaOptions = CobylaOptions(),
)
    options.initial_step > 0.0 || error("initial_step must be positive.")
    options.xtol_rel > 0.0 || error("xtol_rel must be positive.")
    options.maxeval > 0 || error("maxeval must be positive.")
    options.num_starts > 0 || error("num_starts must be positive.")
    options.progress_frequency > 0 || error("progress_frequency must be positive.")

    history = AircraftOptimizationEvaluation[]
    progress_rows = Vector{Vector{Pair{String, Any}}}()
    return_codes = Any[]
    best_candidate = nothing
    best_candidate_model = nothing
    best_feasible = nothing
    best_feasible_model = nothing

    for start in multistart_points(problem, options)
        start_result = run_nlopt_start(problem, start, options)
        push!(return_codes, start_result.return_code)
        append!(history, start_result.history)
        append!(progress_rows, start_result.progress_rows)

        if is_better_evaluation(start_result.best_candidate, best_candidate)
            best_candidate = start_result.best_candidate
            best_candidate_model = start_result.best_candidate_model
        end

        if !isnothing(start_result.best_feasible) &&
           is_better_evaluation(start_result.best_feasible, best_feasible)
            best_feasible = start_result.best_feasible
            best_feasible_model = start_result.best_feasible_model
        end
    end

    optimum = isnothing(best_feasible) ? best_candidate : best_feasible
    optimum_model = isnothing(best_feasible) ? best_candidate_model : best_feasible_model

    return AircraftOptimizationResult(
        isnothing(options.local_algorithm) ?
            options.algorithm :
            Symbol("$(options.algorithm)_$(options.local_algorithm)"),
        return_codes,
        optimum,
        optimum_model,
        best_feasible,
        best_feasible_model,
        history,
        progress_rows,
        length(history),
    )
end

"""
    plot_optimization_history(result, problem; use_scaled=false, columns=3, show_best_objective=true)

Plots the objective and each design variable against evaluation number.
"""
function plot_optimization_history(
    result::AircraftOptimizationResult,
    problem::AircraftOptimizationProblem;
    use_scaled::Bool = false,
    columns::Int = 3,
    show_best_objective::Bool = true,
    show_variable_bounds::Bool = true,
)
    isempty(result.history) && error("No optimization history is available to plot.")
    columns > 0 || error("columns must be positive.")

    eval_index = collect(1:length(result.history))
    objective_history = [evaluation.objective for evaluation in result.history]
    best_objective = accumulate(min, objective_history)

    panels = Any[]
    objective_plot = plot(
        eval_index,
        objective_history;
        xlabel="Evaluation",
        ylabel="Objective",
        title="Objective",
        label="objective",
        linewidth=2,
    )

    if show_best_objective
        plot!(
            objective_plot,
            eval_index,
            best_objective;
            label="best so far",
            linestyle=:dash,
            linewidth=2,
        )
    end

    push!(panels, objective_plot)

    for (ivar, var) in enumerate(problem.design_variables)
        values = if use_scaled
            [evaluation.x_scaled[ivar] for evaluation in result.history]
        else
            [
                design_variable_display_value(var.name, evaluation.x_physical[ivar])
                for evaluation in result.history
            ]
        end

        variable_plot = plot(
            eval_index,
            values;
            xlabel="Evaluation",
            ylabel=use_scaled ? "Scaled Value" : design_variable_display_ylabel(var.name),
            title=use_scaled ? var.name : design_variable_display_title(var.name),
            label=use_scaled ? var.name : design_variable_display_name(var.name),
            linewidth=2,
            legend=false,
        )

        if !use_scaled && show_variable_bounds
            hline!(
                variable_plot,
                [
                    design_variable_display_value(var.name, var.lower),
                    design_variable_display_value(var.name, var.upper),
                ];
                linestyle=:dot,
                color=:gray,
                label="",
            )
        end

        push!(panels, variable_plot)
    end

    nrows = cld(length(panels), columns)
    return plot(
        panels...;
        layout=(nrows, columns),
        size=(450 * columns, 240 * nrows),
    )
end

"""
    weighted_pfei_objective(weights; missions=1:length(weights))

Returns an objective function that forms a weighted sum of mission PFEI values.
"""
function weighted_pfei_objective(
    weights::AbstractVector{<:Real};
    missions = collect(1:length(weights)),
)
    weights_f = Float64.(collect(weights))
    missions_i = Int.(collect(missions))
    length(weights_f) == length(missions_i) || error("weights and missions must have the same length.")

    return ac -> begin
        total = 0.0
        for (weight, mission) in zip(weights_f, missions_i)
            total += weight * ac.parm[imPFEI, mission]
        end
        return total
    end
end

"""
    build_weighted_pfei_problem(base_model, design_variables; ...)

Convenience constructor for weighted multi-mission PFEI optimization.
"""
function build_weighted_pfei_problem(
    base_model::aircraft,
    design_variables::Vector{OptimizationDesignVariable};
    weights::AbstractVector{<:Real},
    missions = collect(1:length(weights)),
    constraints::Vector{OptimizationConstraint} = OptimizationConstraint[],
    mission_modes::Dict{Int, Symbol} = Dict{Int, Symbol}(),
    design_iter::Int = 50,
    offdesign_iter::Int = 50,
    initializes_engine::Bool = true,
    printiter::Bool = false,
    printTO::Bool = false,
    offdesign_missions::Vector{Int} = Int[],
    failure_objective::Real = 1e12,
    failure_constraint::Real = 1e6,
)
    missions_i = Int.(collect(missions))
    nmissions = size(base_model.parm, 2)
    length(weights) == length(missions_i) || error("weights and missions must have the same length.")
    all(1 <= mission <= nmissions for mission in missions_i) || error("missions contains an out-of-bounds mission index.")

    run_settings = MissionOptimizationSettings(
        design_iter=design_iter,
        offdesign_iter=offdesign_iter,
        printiter=printiter,
        printTO=printTO,
        initializes_engine=initializes_engine,
        mission_modes=mission_modes,
        offdesign_missions=offdesign_missions,
    )

    return AircraftOptimizationProblem(
        base_model,
        design_variables,
        weighted_pfei_objective(weights; missions=missions_i);
        constraints=constraints,
        run_settings=run_settings,
        failure_objective=failure_objective,
        failure_constraint=failure_constraint,
    )
end

"""
    balanced_field_length_limit(range_m)

Implements the range-based BFL rule in meters.
"""
function balanced_field_length_limit(range_m::Real)
    range_nmi = Float64(range_m) / nmi_to_m
    if range_nmi < 2500.0
        return 2.2e3
    elseif range_nmi <= 3500.0
        return 2.4e3
    else
        return 3.4e3
    end
end

function fuel_capacity_constraints(; missions = 1:1, tolerance::Real = 1e-3)
    constraints = OptimizationConstraint[]
    for mission in missions
        push!(
            constraints,
            OptimizationConstraint(
                "fuel_capacity_m$(mission)",
                ac -> ac.parm[imWfuel, mission] / ac.parg[igWfmax] - 1.0;
                tolerance=tolerance,
            ),
        )
    end
    return constraints
end

function balanced_field_length_constraints(
    ;
    missions = 1:1,
    tolerance::Real = 1e-3,
    limit_fn::Function = balanced_field_length_limit,
)
    constraints = OptimizationConstraint[]
    for mission in missions
        push!(
            constraints,
            OptimizationConstraint(
                "balanced_field_length_m$(mission)",
                ac -> ac.parm[imlBF, mission] / limit_fn(ac.parm[imRange, mission]) - 1.0;
                tolerance=tolerance,
            ),
        )
    end
    return constraints
end

function wing_span_constraint(; tolerance::Real = 1e-3)
    return OptimizationConstraint(
        "wing_span",
        ac -> ac.wing.layout.span / ac.wing.layout.max_span - 1.0;
        tolerance=tolerance,
    )
end

function top_of_climb_gradient_constraint(; mission::Int = 1, tolerance::Real = 1e-3)
    return OptimizationConstraint(
        "top_of_climb_gradient_m$(mission)",
        ac -> 1.0 - ac.para[iagamV, ipclimbn, mission] / ac.parg[iggtocmin];
        tolerance=tolerance,
    )
end

function tt3_constraint(; mission::Int = 1, limit::Real = 900.0, tolerance::Real = 1e-3)
    limit_f = Float64(limit)
    return OptimizationConstraint(
        "tt3_m$(mission)",
        ac -> maximum(ac.pare[ieTt3, :, mission]) / limit_f - 1.0;
        tolerance=tolerance,
    )
end

function vane_temperature_constraint(; mission::Int = 1, tolerance::Real = 1e-3)
    return OptimizationConstraint(
        "vane_temperature_m$(mission)",
        ac -> maximum(ac.pare[ieTmet1, :, mission]) / ac.parg[igTmetal] - 1.0;
        tolerance=tolerance,
    )
end

function takeoff_weight_constraint(; mission::Int = 1, tolerance::Real = 1e-3)
    return OptimizationConstraint(
        "takeoff_weight_m$(mission)",
        ac -> ac.parm[imWTO, mission] / ac.parg[igWMTO] - 1.0;
        tolerance=tolerance,
    )
end

function fan_diameter_constraint(; limit::Real = 2.0, tolerance::Real = 1e-3)
    limit_f = Float64(limit)
    return OptimizationConstraint(
        "fan_diameter",
        ac -> ac.parg[igdfan] / limit_f - 1.0;
        tolerance=tolerance,
    )
end

function fan_pr_less_than_lpc_pr_constraints(; missions = 1:1, tolerance::Real = 1e-3)
    constraints = OptimizationConstraint[]
    for mission in missions
        push!(
            constraints,
            OptimizationConstraint(
                "fan_pr_le_lpc_pr_m$(mission)",
                ac -> ac.pare[iepif, ipcruise1, mission] / ac.pare[iepilc, ipcruise1, mission] - 1.0;
                tolerance=tolerance,
            ),
        )
    end
    return constraints
end

function standard_weighted_pfei_constraints(
    ;
    missions = 1:1,
    gtoc_mission::Int = 1,
    tt3_mission::Int = 1,
    vane_mission::Int = 1,
    wto_mission::Int = 1,
    tolerance::Real = 1e-3,
    tt3_limit::Real = 900.0,
    fan_diameter_limit::Real = 2.0,
)
    constraints = OptimizationConstraint[]
    append!(constraints, fuel_capacity_constraints(; missions=missions, tolerance=tolerance))
    append!(constraints, balanced_field_length_constraints(; missions=missions, tolerance=tolerance))
    push!(constraints, wing_span_constraint(; tolerance=tolerance))
    push!(constraints, top_of_climb_gradient_constraint(; mission=gtoc_mission, tolerance=tolerance))
    push!(constraints, tt3_constraint(; mission=tt3_mission, limit=tt3_limit, tolerance=tolerance))
    push!(constraints, vane_temperature_constraint(; mission=vane_mission, tolerance=tolerance))
    push!(constraints, takeoff_weight_constraint(; mission=wto_mission, tolerance=tolerance))
    push!(constraints, fan_diameter_constraint(; limit=fan_diameter_limit, tolerance=tolerance))
    append!(constraints, fan_pr_less_than_lpc_pr_constraints(; missions=missions, tolerance=tolerance))
    return constraints
end

function default_initial(lower::Real, upper::Real, initial::Union{Nothing, Real})
    return isnothing(initial) ? 0.5 * (Float64(lower) + Float64(upper)) : Float64(initial)
end

function aspect_ratio_variable(; lower::Real = 7.0, upper::Real = 12.0, initial::Union{Nothing, Real} = nothing)
    return OptimizationDesignVariable(
        "wing_aspect_ratio",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.wing.layout.AR = value),
    )
end

function wing_sweep_variable(; lower::Real = 10.0, upper::Real = 40.0, initial::Union{Nothing, Real} = nothing)
    return OptimizationDesignVariable(
        "wing_sweep",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.wing.layout.sweep = value),
    )
end

function inboard_taper_ratio_variable(; lower::Real = 0.1, upper::Real = 1.0, initial::Union{Nothing, Real} = nothing)
    return OptimizationDesignVariable(
        "inboard_taper_ratio",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.wing.inboard.λ = value),
    )
end

function outboard_taper_ratio_variable(; lower::Real = 0.1, upper::Real = 1.0, initial::Union{Nothing, Real} = nothing)
    return OptimizationDesignVariable(
        "outboard_taper_ratio",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.wing.outboard.λ = value),
    )
end

function root_thickness_to_chord_variable(
    ;
    lower::Real = 0.1,
    upper::Real = 0.15,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "root_thickness_to_chord",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.wing.inboard.cross_section.thickness_to_chord = value),
    )
end

function spanbreak_thickness_to_chord_variable(
    ;
    lower::Real = 0.1,
    upper::Real = 0.15,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "spanbreak_thickness_to_chord",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.wing.outboard.cross_section.thickness_to_chord = value),
    )
end

function mission_cruise_altitude_variable(
    mission::Int;
    lower::Real = 20_000.0 * ft_to_m,
    upper::Real = 60_000.0 * ft_to_m,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "cruise_altitude_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.para[iaalt, ipcruise1, mission] = value),
    )
end

function mission_cruise_cl_variable(
    mission::Int;
    lower::Real = 0.45,
    upper::Real = 0.65,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "cruise_cl_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.para[iaCL, ipclimb1+1:ipdescentn-1, mission] .= value),
    )
end

function mission_rcls_variable(
    mission::Int;
    lower::Real = 0.1,
    upper::Real = 1.4,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "rcls_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.para[iarcls, ipclimb1+1:ipdescentn-1, mission] .= value),
    )
end

function mission_rclt_variable(
    mission::Int;
    lower::Real = 0.1,
    upper::Real = 1.4,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "rclt_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.para[iarclt, ipclimb1+1:ipdescentn-1, mission] .= value),
    )
end

function shared_mission_rcls_variable(
    ;
    lower::Real = 0.1,
    upper::Real = 1.4,
    initial::Union{Nothing, Real} = nothing,
    missions::AbstractVector{<:Integer} = [1, 2, 3],
    name::AbstractString = "rcls_shared",
)
    missions_i = Int.(collect(missions))
    isempty(missions_i) && error("missions must contain at least one mission index.")
    return OptimizationDesignVariable(
        name,
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.para[iarcls, ipclimb1+1:ipdescentn-1, missions_i] .= value),
    )
end

function shared_mission_rclt_variable(
    ;
    lower::Real = 0.1,
    upper::Real = 1.4,
    initial::Union{Nothing, Real} = nothing,
    missions::AbstractVector{<:Integer} = [1, 2, 3],
    name::AbstractString = "rclt_shared",
)
    missions_i = Int.(collect(missions))
    isempty(missions_i) && error("missions must contain at least one mission index.")
    return OptimizationDesignVariable(
        name,
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.para[iarclt, ipclimb1+1:ipdescentn-1, missions_i] .= value),
    )
end

function shared_cruise_turbine_temperature_variable(
    ;
    lower::Real = 700.0,
    upper::Real = 2000.0,
    initial::Union{Nothing, Real} = nothing,
    missions::Union{Nothing, AbstractVector{Int}} = nothing,
)
    setter! = if isnothing(missions)
        (ac, value) -> (ac.pare[ieTt4, ipcruise1:ipcruise2, :] .= value)
    else
        missions_i = Int.(collect(missions))
        (ac, value) -> (ac.pare[ieTt4, ipcruise1:ipcruise2, missions_i] .= value)
    end

    return OptimizationDesignVariable(
        "shared_cruise_turbine_temperature",
        lower,
        upper,
        default_initial(lower, upper, initial),
        setter!,
    )
end

function shared_takeoff_turbine_temperature_variable(
    ;
    lower::Real = 1100.0,
    upper::Real = 2000.0,
    initial::Union{Nothing, Real} = nothing,
    missions::Union{Nothing, AbstractVector{Int}} = nothing,
)
    setter! = if isnothing(missions)
        (ac, value) -> (ac.pare[ieTt4, ipstatic:iptakeoff, :] .= value)
    else
        missions_i = Int.(collect(missions))
        (ac, value) -> (ac.pare[ieTt4, ipstatic:iptakeoff, missions_i] .= value)
    end

    return OptimizationDesignVariable(
        "shared_takeoff_turbine_temperature",
        lower,
        upper,
        default_initial(lower, upper, initial),
        setter!,
    )
end

function shared_fan_pr_variable(
    ;
    lower::Real = 1.0,
    upper::Real = 3.0,
    initial::Union{Nothing, Real} = nothing,
    missions::Union{Nothing, AbstractVector{Int}} = nothing,
)
    setter! = if isnothing(missions)
        (ac, value) -> (ac.pare[iepif, ipcruise1, :] .= value)
    else
        missions_i = Int.(collect(missions))
        (ac, value) -> (ac.pare[iepif, ipcruise1, missions_i] .= value)
    end

    return OptimizationDesignVariable(
        "shared_fan_pr",
        lower,
        upper,
        default_initial(lower, upper, initial),
        setter!,
    )
end

function shared_lpc_pr_variable(
    ;
    lower::Real = 1.0,
    upper::Real = 50.0,
    initial::Union{Nothing, Real} = nothing,
    missions::Union{Nothing, AbstractVector{Int}} = nothing,
)
    setter! = if isnothing(missions)
        (ac, value) -> (ac.pare[iepilc, ipcruise1, :] .= value)
    else
        missions_i = Int.(collect(missions))
        (ac, value) -> (ac.pare[iepilc, ipcruise1, missions_i] .= value)
    end

    return OptimizationDesignVariable(
        "shared_lpc_pr",
        lower,
        upper,
        default_initial(lower, upper, initial),
        setter!,
    )
end

function shared_hpc_pr_variable(
    ;
    lower::Real = 1.0,
    upper::Real = 50.0,
    initial::Union{Nothing, Real} = nothing,
    missions::Union{Nothing, AbstractVector{Int}} = nothing,
)
    setter! = if isnothing(missions)
        (ac, value) -> (ac.pare[iepihc, ipcruise1, :] .= value)
    else
        missions_i = Int.(collect(missions))
        (ac, value) -> (ac.pare[iepihc, ipcruise1, missions_i] .= value)
    end

    return OptimizationDesignVariable(
        "shared_hpc_pr",
        lower,
        upper,
        default_initial(lower, upper, initial),
        setter!,
    )
end

function shared_bpr_variable(
    ;
    lower::Real = 1.0,
    upper::Real = 12.0,
    initial::Union{Nothing, Real} = nothing,
    missions::Union{Nothing, AbstractVector{Int}} = nothing,
)
    setter! = if isnothing(missions)
        (ac, value) -> (ac.pare[ieBPR, ipcruise1, :] .= value)
    else
        missions_i = Int.(collect(missions))
        (ac, value) -> (ac.pare[ieBPR, ipcruise1, missions_i] .= value)
    end

    return OptimizationDesignVariable(
        "shared_bpr",
        lower,
        upper,
        default_initial(lower, upper, initial),
        setter!,
    )
end

function deviation_height_variable(
    mission::Int = 3;
    lower::Real = -2133.6,
    upper::Real = 2133.6,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "deviation_height_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.parm[imDeviationDH, mission] = value),
    )
end

function deviation_length_variable(
    mission::Int = 3;
    lower::Real = 0.0,
    upper::Real,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "deviation_length_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.parm[imDeviationDL, mission] = value),
    )
end

function deviation_start_variable(
    mission::Int = 3;
    lower::Real = 0.0,
    upper::Real,
    initial::Union{Nothing, Real} = nothing,
)
    return OptimizationDesignVariable(
        "deviation_start_m$(mission)",
        lower,
        upper,
        default_initial(lower, upper, initial),
        (ac, value) -> (ac.parm[imDeviationStart, mission] = value),
    )
end
