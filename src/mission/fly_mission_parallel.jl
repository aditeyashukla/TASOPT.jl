"""
    fly_mission_parallel!(ac, missions; kwargs...)

Fly independent off-design missions in parallel after the aircraft has already
been sized. Each mission is evaluated on a deep copy of `ac`, and only that
mission's results are merged back into `ac`.

The parallel work uses Julia `Distributed` worker processes. If no worker
processes are available, workers are added by default.

# Keywords
- `itermax=35`
- `initializes_engine=true`
- `opt_prescribed_cruise_parameter="altitude"`
- `printTO=false`
- `nworkers=min(length(missions), Sys.CPU_THREADS)`
- `add_workers=true`
- `cleanup_workers=false`
- `exeflags=nothing`
- `collect_deviation_diagnostics=false`
"""
function fly_mission_parallel!(
    ac::aircraft,
    missions;
    itermax::Integer = 35,
    initializes_engine::Bool = true,
    opt_prescribed_cruise_parameter = "altitude",
    printTO::Bool = false,
    wrlx1::Real = 0.5,
    wrlx2::Real = 0.9,
    wrlx3::Real = 0.5,
    maximum_deviation_surrogates::Integer = typemax(Int),
    allow_deviation_t4_violation::Bool = false,
    nworkers::Integer = _default_parallel_worker_count(missions),
    add_workers::Bool = true,
    cleanup_workers::Bool = false,
    exeflags = nothing,
    collect_deviation_diagnostics::Bool = false,
)
    mission_ids = _normalize_parallel_missions(ac, missions)

    if isempty(ac.is_sized) || !ac.is_sized[1]
        error("Aircraft not sized. Please size the aircraft before running off-design missions.")
    end

    mission_kwargs = (
        itermax = Int(itermax),
        initializes_engine = initializes_engine,
        opt_prescribed_cruise_parameter = opt_prescribed_cruise_parameter,
        printTO = printTO,
        wrlx1 = Float64(wrlx1),
        wrlx2 = Float64(wrlx2),
        wrlx3 = Float64(wrlx3),
        maximum_deviation_surrogates = Int(maximum_deviation_surrogates),
        allow_deviation_t4_violation = allow_deviation_t4_violation,
    )

    worker_count = max(1, min(Int(nworkers), length(mission_ids)))
    pids, added_pids = _ensure_parallel_mission_workers!(
        worker_count;
        add_workers = add_workers,
        exeflags = exeflags,
    )

    try
        _load_tasopt_on_workers!(pids)

        pool = Distributed.CachingPool(pids)
        if collect_deviation_diagnostics
            for imission in mission_ids
                _clear_deviation_diagnostic!(ac, imission)
            end
        end
        results = Distributed.pmap(pool, mission_ids) do imission
            _fly_mission_worker_result(
                ac, imission, mission_kwargs, collect_deviation_diagnostics,
            )
        end

        for result in results
            _merge_parallel_mission_result!(ac, result)
        end

        return ac
    finally
        if cleanup_workers && !isempty(added_pids)
            Distributed.rmprocs(added_pids...)
            setdiff!(_tasopt_loaded_parallel_workers, added_pids)
        end
    end
end

# Backward-compatible alias for the spelling used in existing notes/examples.
fly_mission_paralell!(args...; kwargs...) = fly_mission_parallel!(args...; kwargs...)

"""
    size_aircraft_and_fly_missions_parallel!(ac, missions; kwargs...)

Size `ac` and then fly the requested off-design missions in parallel.

Sizing keywords are prefixed with `size_`; mission and worker keywords match
`fly_mission_parallel!`.
"""
function size_aircraft_and_fly_missions_parallel!(
    ac::aircraft,
    missions;
    size_iter::Integer = 35,
    initwgt::Bool = false,
    Ldebug::Bool = false,
    printiter::Bool = true,
    itermax::Integer = 35,
    initializes_engine::Bool = true,
    opt_prescribed_cruise_parameter = "altitude",
    printTO::Bool = false,
    wrlx1::Real = 0.5,
    wrlx2::Real = 0.9,
    wrlx3::Real = 0.5,
    maximum_deviation_surrogates::Integer = typemax(Int),
    allow_deviation_t4_violation::Bool = false,
    nworkers::Integer = _default_parallel_worker_count(missions),
    add_workers::Bool = true,
    cleanup_workers::Bool = false,
    exeflags = nothing,
    collect_deviation_diagnostics::Bool = false,
)
    size_aircraft!(
        ac;
        iter = Int(size_iter),
        initwgt = initwgt,
        Ldebug = Ldebug,
        printiter = printiter,
    )

    return fly_mission_parallel!(
        ac,
        missions;
        itermax = itermax,
        initializes_engine = initializes_engine,
        opt_prescribed_cruise_parameter = opt_prescribed_cruise_parameter,
        printTO = printTO,
        wrlx1 = wrlx1,
        wrlx2 = wrlx2,
        wrlx3 = wrlx3,
        maximum_deviation_surrogates = maximum_deviation_surrogates,
        allow_deviation_t4_violation = allow_deviation_t4_violation,
        nworkers = nworkers,
        add_workers = add_workers,
        cleanup_workers = cleanup_workers,
        exeflags = exeflags,
        collect_deviation_diagnostics = collect_deviation_diagnostics,
    )
end

function _normalize_parallel_missions(ac::aircraft, missions)
    mission_ids = _parallel_mission_ids(missions)
    isempty(mission_ids) && error("At least one mission index is required.")
    length(unique(mission_ids)) == length(mission_ids) ||
        error("Mission indices must be unique.")

    nmissions = size(ac.parm, 2)
    for imission in mission_ids
        1 <= imission <= nmissions ||
            error("Mission index $(imission) is outside the available mission range 1:$(nmissions).")
    end

    return mission_ids
end

_parallel_mission_ids(mission::Integer) = [Int(mission)]
_parallel_mission_ids(missions) = Int.(collect(missions))

function _default_parallel_worker_count(missions)
    return min(length(_parallel_mission_ids(missions)), Sys.CPU_THREADS)
end

function _ensure_parallel_mission_workers!(
    requested_workers::Integer;
    add_workers::Bool,
    exeflags,
)
    requested_workers = max(1, Int(requested_workers))
    remote_workers = Distributed.nprocs() == 1 ? Int[] : Distributed.workers()
    added_pids = Int[]

    needed = requested_workers - length(remote_workers)
    if needed > 0 && add_workers
        added_pids = if isnothing(exeflags)
            Distributed.addprocs(needed)
        else
            Distributed.addprocs(needed; exeflags = exeflags)
        end
        remote_workers = Distributed.workers()
    end

    pids = isempty(remote_workers) ? [Distributed.myid()] : remote_workers
    return pids, added_pids
end

function _load_tasopt_on_workers!(pids)
    remote_pids = setdiff(collect(pids), [Distributed.myid()])
    isempty(remote_pids) && return nothing

    unloaded_pids = setdiff(remote_pids, collect(_tasopt_loaded_parallel_workers))
    isempty(unloaded_pids) && return nothing

    Distributed.remotecall_eval(Main, unloaded_pids, :(using TASOPT))
    union!(_tasopt_loaded_parallel_workers, unloaded_pids)
    return nothing
end

const _tasopt_loaded_parallel_workers = Set{Int}()

function _fly_mission_worker_result(
    ac_sized,
    imission::Integer,
    mission_kwargs::NamedTuple,
    collect_deviation_diagnostics::Bool = false,
)
    ac_work = deepcopy(ac_sized)
    if !collect_deviation_diagnostics
        fly_mission!(ac_work, Int(imission); mission_kwargs...)
        return _parallel_mission_result(ac_work, Int(imission))
    end

    diagnostics_were_enabled = _deviation_diagnostics_enabled[]
    set_deviation_diagnostics!(true)
    clear_deviation_diagnostics!()
    try
        fly_mission!(ac_work, Int(imission); mission_kwargs...)
        diagnostic = deviation_diagnostics(ac_work, Int(imission))
        return merge(
            _parallel_mission_result(ac_work, Int(imission)),
            (; deviation_diagnostic = diagnostic),
        )
    finally
        clear_deviation_diagnostics!()
        set_deviation_diagnostics!(diagnostics_were_enabled)
    end
end

function _parallel_mission_result(ac::aircraft, imission::Int)
    return (
        imission = imission,
        parm = copy(ac.parm[:, imission]),
        para = copy(ac.para[:, :, imission]),
        pare = copy(ac.pare[:, :, imission]),
        engine_state = _extract_engine_mission_state(ac.engine, imission, size(ac.parm, 2)),
    )
end

function _merge_parallel_mission_result!(ac::aircraft, result)
    imission = result.imission
    ac.parm[:, imission] .= result.parm
    ac.para[:, :, imission] .= result.para
    ac.pare[:, :, imission] .= result.pare
    _merge_engine_mission_state!(ac.engine, result.engine_state, imission, size(ac.parm, 2))
    if hasproperty(result, :deviation_diagnostic) && !isnothing(result.deviation_diagnostic)
        mission_data = get!(_deviation_diagnostics, ac, Dict{Int, NamedTuple}())
        mission_data[imission] = result.deviation_diagnostic
    end
    return ac
end

function _extract_engine_mission_state(src_engine, imission::Int, nmissions::Int)
    return (
        heat_exchangers = [
            _extract_heat_exchanger_mission_state(src_hx, imission)
            for src_hx in src_engine.heat_exchangers
        ],
        data = _extract_mission_indexed_data_fields(src_engine.data, imission, nmissions),
    )
end

function _merge_engine_mission_state!(dest_engine, src_state, imission::Int, nmissions::Int)
    for (dest_hx, src_hx_state) in zip(dest_engine.heat_exchangers, src_state.heat_exchangers)
        _copy_heat_exchanger_mission_state!(dest_hx, src_hx_state, imission)
    end

    _copy_mission_indexed_data_fields!(dest_engine.data, src_state.data, imission, nmissions)
    return dest_engine
end

function _extract_heat_exchanger_mission_state(src_hx, imission::Int)
    src_mission = src_hx.HXgas_mission

    if ndims(src_mission) != 2 || size(src_mission, 2) < imission
        return nothing
    end

    npoints = size(src_mission, 1)
    assigned = falses(npoints)
    values = Vector{Any}(undef, npoints)
    for ip in 1:npoints
        if isassigned(src_mission, ip, imission)
            values[ip] = deepcopy(src_mission[ip, imission])
            assigned[ip] = true
        end
    end

    return (assigned = assigned, values = values)
end

function _copy_heat_exchanger_mission_state!(dest_hx, src_hx_state, imission::Int)
    src_hx_state === nothing && return dest_hx

    dest_mission = dest_hx.HXgas_mission
    if ndims(dest_mission) != 2 || size(dest_mission, 2) < imission
        return dest_hx
    end

    npoints = min(size(dest_mission, 1), length(src_hx_state.assigned))
    for ip in 1:npoints
        if src_hx_state.assigned[ip]
            dest_mission[ip, imission] = src_hx_state.values[ip]
        end
    end

    return dest_hx
end

function _extract_mission_indexed_data_fields(src_data, imission::Int, nmissions::Int)
    Base.ismutabletype(typeof(src_data)) || return NamedTuple[]

    fields = NamedTuple[]
    for field in fieldnames(typeof(src_data))
        src_value = getfield(src_data, field)
        if src_value isa AbstractMatrix &&
                size(src_value, 2) == nmissions &&
                size(src_value, 2) >= imission
            push!(fields, (field = field, value = copy(src_value[:, imission])))
        end
    end

    return fields
end

function _copy_mission_indexed_data_fields!(dest_data, src_data_fields, imission::Int, nmissions::Int)
    Base.ismutabletype(typeof(dest_data)) || return dest_data

    for src_field in src_data_fields
        hasfield(typeof(dest_data), src_field.field) || continue

        dest_value = getfield(dest_data, src_field.field)
        if dest_value isa AbstractMatrix &&
                size(dest_value, 2) == nmissions &&
                size(dest_value, 2) >= imission &&
                size(dest_value, 1) == length(src_field.value)
            dest_value[:, imission] .= src_field.value
        end
    end

    return dest_data
end
