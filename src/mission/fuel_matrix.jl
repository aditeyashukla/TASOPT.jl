"""
    avoidance_fuel_matrix(ac, mission; cruise_altitudes_ft, deviation_heights_ft,
                          deviation_lengths_nmi, start_nmi=0, itermax=600,
                          cell_callback=nothing, mission_weight_rtol=1e-8, integration_rtol=1e-7, integration_max_step=500*nmi_to_m)

Replay a sized aircraft without resizing or modifying `ac`. Returns `fuel_kg`
(trip fuel excluding reserve), `loaded_fuel_kg`, `trip_burn_kg`, `feasible`, and diagnostics
on axes (cruise altitude, signed DH, DL). Failed cells are NaN, never zero.
Either zero deviation resource is the canonical no-maneuver mission.

Each cell starts from an independent copy of the input aircraft. The maneuver
starts at `start_nmi` after start of cruise; this is a shared-placement fuel
approximation when combined with best-start weather matrices. Serial use only:
TASOPT's deviation telemetry uses process-global state.
"""
function avoidance_fuel_matrix(ac, mission::Integer;
        cruise_altitudes_ft=[ac.para[iaalt, ipcruise1, mission] / ft_to_m],
        deviation_heights_ft, deviation_lengths_nmi, start_nmi=0., itermax=600,
        cell_callback=nothing, mission_weight_rtol=1e-8, integration_rtol=1e-7, integration_max_step=500*nmi_to_m)
    ac.is_sized[1] || throw(ArgumentError("Aircraft must already be sized"))
    1 < mission <= size(ac.parm, 2) || throw(ArgumentError("Expected an off-design mission"))
    axes = (Float64.(collect(cruise_altitudes_ft)), Float64.(collect(deviation_heights_ft)),
            Float64.(collect(deviation_lengths_nmi)))
    all(a -> !isempty(a) && all(isfinite, a) && length(unique(a)) == length(a), axes) ||
        throw(ArgumentError("Axes must be nonempty, finite and unique"))
    all(>=(0), axes[3]) && isfinite(start_nmi) && start_nmi >= 0 ||
        throw(ArgumentError("Lengths and start must be nonnegative"))
    itermax > 0 || throw(ArgumentError("itermax must be positive"))
    fuel = fill(NaN, length.(axes)...)
    feasible = falses(size(fuel))
    diagnostics = Array{Any}(undef, size(fuel))
    cache = Dict{NTuple{3,Float64},Any}()
    prior_diagnostics = _deviation_diagnostics_enabled[]
    set_deviation_diagnostics!(true)
    try
        for ih in eachindex(axes[1]), j in eachindex(axes[2]), k in eachindex(axes[3])
            h, dh, dl = axes[1][ih], axes[2][j], axes[3][k]
            key = dh == 0 || dl == 0 ? (h, 0., 0.) : (h, dh, dl)
            if !haskey(cache, key)
                trial = deepcopy(ac)
                trial.para[iaalt, ipcruise1, mission] = h * ft_to_m
                trial.parm[imDeviationDH, mission] = key[2] * ft_to_m
                trial.parm[imDeviationDL, mission] = key[3] * nmi_to_m
                trial.parm[imDeviationStart, mission] = start_nmi * nmi_to_m
                telemetry = Ref{Any}(nothing)
                started = time()
                fb = NaN; residual = NaN; iterations = 0
                tmargin = NaN; gmargin = NaN; remaining = NaN
                error_message = ""; ok = false
                try
                    redirect_stdout(devnull) do
                        fly_mission!(trial, mission; itermax, printTO=false,
                            initializes_engine=true, opt_prescribed_cruise_parameter="altitude",
                            allow_deviation_t4_violation=false,
                            mission_weight_rtol, integration_rtol, integration_max_step,
                            iteration_callback=s -> (telemetry[] = s))
                    end
                    t = telemetry[]
                    t === nothing && error("Missing convergence telemetry")
                    residual = t.errw; iterations = t.iteration
                    residual <= 1e-8 && !t.used_deviation_surrogate ||
                        error("Unconverged mission or deviation surrogate")
                    all(trial.pare[ieConvFail, :, mission] .== 0) || error("Engine convergence failure")
                    isapprox(trial.para[iaalt, ipcruise1, mission], h*ft_to_m; atol=1e-6, rtol=0) ||
                        error("Cruise altitude changed")
                    if key[2] != 0 && key[3] != 0
                        d = deviation_diagnostics(trial, mission)
                        d === nothing && error("Missing deviation diagnostics")
                        isapprox(first(d.local_range_m), start_nmi*nmi_to_m; atol=1e-6, rtol=0) ||
                            error("Maneuver start differs")
                        tmargin = d.t4_feasibility_margin; gmargin = d.gamma_feasibility_margin
                        tmargin <= 0 && gmargin <= 1e-7 || error("Infeasible deviation temperature or climb angle")
                        remaining = (trial.para[iaRange, ipcruisen, mission] -
                                     trial.para[iaRange, ipdeviation4, mission]) / nmi_to_m
                        remaining >= -1e-7 || error("Deviation extends beyond cruise")
                    end
                    fb = trial.parm[imWfuel, mission] / gee
                    isfinite(fb) && fb > 0 || error("Invalid mission fuel")
                    fb / (trial.parg[igWfmax]/gee) - 1 <= 1e-6 || error("Fuel capacity exceeded")
                    ok = true
                catch err
                    err isa Union{InterruptException,MethodError,UndefVarError,BoundsError,TypeError} && rethrow()
                    error_message = replace(sprint(showerror, err), '\n' => ' ')
                    fb = NaN
                finally
                    delete!(_deviation_diagnostics, trial)
                end
                cache[key] = (; success=ok, fuel_kg=fb/(1+ac.parg[igfreserve]), loaded_fuel_kg=fb,
                    trip_burn_kg=fb/(1+ac.parg[igfreserve]), residual, iterations,
                    temperature_margin_K=tmargin, gamma_margin=gmargin,
                    remaining_cruise_nmi=remaining, seconds=time()-started, error=error_message)
            end
            r = cache[key]
            fuel[ih,j,k] = r.fuel_kg; feasible[ih,j,k] = r.success; diagnostics[ih,j,k] = r
            cell_callback === nothing || cell_callback((; h_ft=h, dH_ft=dh, dL_nmi=dl, r...))
        end
    finally
        set_deviation_diagnostics!(prior_diagnostics)
    end
    (; fuel_kg=fuel, trip_burn_kg=copy(fuel), loaded_fuel_kg=fuel.*(1+ac.parg[igfreserve]), feasible, diagnostics,
       cruise_altitudes_ft=axes[1], deviation_heights_ft=axes[2], deviation_lengths_nmi=axes[3])
end
