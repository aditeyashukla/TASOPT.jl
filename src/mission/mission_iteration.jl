"""
    _mission_iteration!(ac, imission, Ldebug; calculate_cruise = false)

Runs aircraft through mission, calculating fuel burn
and other mission variables. Formerly, `mission!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft data storage object
    - `imission::Int64`: mission index
    - `Ldebug::Bool`: debugging flag

NOTE: 
This routine assumes that estimates of the climb-leg flight path 
gamma angles are passed in via para[iagamV,ipclimb1:ipclimbn].
These appear as cos(gamma) factors in the climb equations,
and can be passed in as zero with only a minor error.
They are updated and returned in the same para[iagamV,ip] array.

"""
const _deviation_diagnostics_enabled = Ref(false)
const _deviation_diagnostics = IdDict{Any, Dict{Int, NamedTuple}}()

function set_deviation_diagnostics!(enabled::Bool = true)
      _deviation_diagnostics_enabled[] = enabled
      return enabled
end

function clear_deviation_diagnostics!()
      empty!(_deviation_diagnostics)
      return nothing
end

function deviation_diagnostics(ac, imission::Integer)
      mission_data = get(_deviation_diagnostics, ac, nothing)
      mission_data === nothing && return nothing
      return get(mission_data, Int(imission), nothing)
end

function _clear_deviation_diagnostic!(ac, imission::Integer)
      _deviation_diagnostics_enabled[] || return nothing
      mission_data = get(_deviation_diagnostics, ac, nothing)
      mission_data === nothing && return nothing
      delete!(mission_data, Int(imission))
      return nothing
end

function _store_deviation_diagnostic!(
      ac,
      imission::Integer,
      range_start,
      dev_s,
      dev_alt,
      dev_time,
      dev_fracW,
      dev_gamma,
      dev_mode,
      gamma_geometry_min,
      gamma_thrust_limit,
      t4_geometry_required,
      t4_avoidance_limit,
      t4_limit_tolerance,
      dev_para,
      dev_pare,
      dev_FFC,
      dev_Vgi,
      dev_FoW,
)
      _deviation_diagnostics_enabled[] || return nothing
      mission_data = get!(_deviation_diagnostics, ac, Dict{Int, NamedTuple}())
      mission_data[Int(imission)] = (
            label = ["D$(i)" for i in eachindex(dev_s)],
            local_range_m = copy(dev_s),
            range_m = range_start .+ copy(dev_s),
            altitude_m = copy(dev_alt),
            time_s = copy(dev_time),
            fracW = copy(dev_fracW),
            gamma = copy(dev_gamma),
            mode = String.(dev_mode),
            gamma_geometry_min = gamma_geometry_min,
            gamma_thrust_limit = gamma_thrust_limit,
            gamma_feasibility_margin = gamma_geometry_min - gamma_thrust_limit,
            t4_geometry_required = t4_geometry_required,
            t4_avoidance_limit = t4_avoidance_limit,
            t4_limit_tolerance = t4_limit_tolerance,
            t4_feasibility_margin = maximum(dev_pare[ieTt4,findall(==(:climb),dev_mode)]) -
                (t4_avoidance_limit + t4_limit_tolerance),
            para = copy(dev_para),
            pare = copy(dev_pare),
            FFC = copy(dev_FFC),
            Vgi = copy(dev_Vgi),
            FoW = copy(dev_FoW),
      )
      return nothing
end

function _mission_iteration!(
      ac,
      imission,
      Ldebug;
      calculate_cruise = false,
      allow_deviation_t4_violation = false,
      integration_rtol = 1e-7,
      integration_cache = true,
      integration_max_step = 500*nmi_to_m,
)
      #Unpack aircraft
      parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, eng, landing_gear = unpack_ac(ac, imission) 
      _clear_deviation_diagnostic!(ac, imission)

      ifirst = true

      # HACK TODO add the para back
      # iairf
      initializes_engine = true

      # Mission range
      Rangetot = parm[imRange]
      para[iaRange, ipdescentn] = Rangetot

      # MTOW
      WMTO = parg[igWMTO]

      # Payload fraction for this mission
      rpay = parm[imWpay] / parg[igWpay]
      ξpay = 0.0

      # Zero-fuel weight for this mission
      Wzero = WMTO - parg[igWfuel] - parg[igWpay] + parm[imWpay]

      # mission TO weight
      WTO = parm[imWTO]

      # mission TO fuel weight
      WfTO = WTO - Wzero

      # set known operating conditions

      # takeoff altitude conditions
      ΔTatmos = parm[imDeltaTatm]
      ip = ipstatic
      alt_m = para[iaalt, ip]
      atmos_state = atmos(alt_m, ΔTatmos)
      T0 = atmos_state.T
      p0 = atmos_state.p
      ρ0 = atmos_state.ρ
      a0 = atmos_state.a
      μ0 = atmos_state.μ

      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = ρ0
      pare[iemu0, ip] = μ0
      pare[ieM0, ip] = 0.0
      pare[ieu0, ip] = 0.0

      para[iaMach, ip] = 0.0
      para[iaCL, ip] = 0.0
      para[iaReunit, ip] = 0.0
      para[iaWbuoy, ip] = 0.0

      # ip ∈ iprotate:ipclimb1
      pare[iep0, iprotate:ipclimb1] .= p0
      pare[ieT0, iprotate:ipclimb1] .= T0
      pare[iea0, iprotate:ipclimb1] .= a0
      pare[ierho0, iprotate:ipclimb1] .= ρ0
      pare[iemu0, iprotate:ipclimb1] .= μ0
      para[iaWbuoy, iprotate:ipclimb1] .= 0.0


      # Start-of-cruise altitude conditions
      ip = ipcruise1
      Mach = para[iaMach, ip]
      alt_m = para[iaalt, ip]

      atmos_state = atmos(alt_m, ΔTatmos)
      T0 = atmos_state.T
      p0 = atmos_state.p
      rho0 = atmos_state.ρ
      a0 = atmos_state.a
      mu0 = atmos_state.μ
      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = rho0
      pare[iemu0, ip] = mu0
      pare[ieM0, ip] = Mach
      pare[ieu0, ip] = Mach * a0
      para[iaReunit, ip] = Mach * a0 * rho0 / mu0

      # End-of-descent altitude conditions
      ip = ipdescentn
      alt_m = para[iaalt, ip]
      atmos_state = atmos(alt_m, ΔTatmos)
      T0 = atmos_state.T
      p0 = atmos_state.p
      rho0 = atmos_state.ρ
      a0 = atmos_state.a
      mu0 = atmos_state.μ

      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = rho0
      pare[iemu0, ip] = mu0
      para[iaWbuoy, ip] = 0.0

      # interpolate CL over climb points, 
      #  between specified ipclimb1+1, ipclimbn values
      CLa = para[iaCL, ipclimb1+1]
      CLb = para[iaCL, ipcruise1]

      @inbounds for ip = ipclimb1+1:ipclimbn
            frac = float(ip - (ipclimb1 + 1)) /
                   float(ipclimbn - (ipclimb1 + 1))
            para[iaCL, ip] = CLa * (1.0 - frac^2) +
                             CLb * frac^2
      end

      #---- interpolate CL over descent points, 
      #-     between specified ipdescent1, ipdescentn-1 values
      #      CLd = para[iaCL,ipcruisen]
      #      CLe = para[iaCL,ipcruisen]
      CLd = para[iaCL, ipcruisen] * 0.96
      CLe = para[iaCL, ipcruisen] * 0.50

      @inbounds for ip = ipdescent1:ipdescentn-1
            frac = float(ip - ipdescent1) /
                   float((ipdescentn - 1) - ipdescent1)
            fb = 1.0 - frac
            para[iaCL, ip] = CLd * fb^2 +
                             CLe * (1.0 - fb^2)
      end


      #---- estimate takeoff speed and set V,Re over climb and descent
      cosL = cosd(wing.layout.sweep)
      CLTO = para[iaclpmax, iptakeoff] * cosL^2
      # [prash] I think the assumption here is that Wcruise/WTO~1 and
      # just scaling it with rho and CL to get an initial estimate for V
      VTO = pare[ieu0, ipcruise1] *
            sqrt(pare[ierho0, ipcruise1] / pare[ierho0, iptakeoff]) *
            sqrt(para[iaCL, ipcruise1] / CLTO)
      ReTO = VTO * pare[ierho0, iptakeoff] / pare[iemu0, iptakeoff]
      pare[ieu0, iprotate:ipclimb1] .= VTO
      para[iaReunit, iprotate:ipclimb1] .= ReTO
      @inbounds for ip = ipclimb1+1:ipclimbn
            frac = float(ip - ipclimb1) / float(ipclimbn - ipclimb1)
            V = VTO * (1.0 - frac) + pare[ieu0, ipcruise1] * frac
            Re = ReTO * (1.0 - frac) + para[iaReunit, ipcruise1] * frac
            pare[ieu0, ip] = V
            para[iaReunit, ip] = Re
      end
      @inbounds for ip = ipdescent1:ipdescentn
            frac = float(ip - ipdescent1) / float(ipdescentn - ipdescent1)
            V = VTO * frac + pare[ieu0, ipcruisen] * (1.0 - frac)
            Re = ReTO * frac + para[iaReunit, ipcruisen] * (1.0 - frac)
            pare[ieu0, ip] = V
            para[iaReunit, ip] = Re
      end



      #---- takeoff CLmax via section clmax and sweep correction
      ip = iprotate
      clpmax = para[iaclpmax, ip]
      sweep = wing.layout.sweep
      cosL = cosd(sweep)
      CLmax = clpmax * cosL^2

      #---- Vs stall speed (takeoff condition)
      rho0 = pare[ierho0, ip]
      a0 = pare[iea0, ip]
      S = wing.layout.S
      Vstall = sqrt(2.0 * WTO / (rho0 * S * CLmax))
      Mstall = Vstall / a0
      pare[ieu0, ip] = Vstall
      pare[ieM0, ip] = Mstall
      para[iaMach, ip] = Mstall
      para[iaReunit, ip] = Vstall * pare[ierho0, ip] / pare[iemu0, ip]

      #---- V2 speed per FAR-25  (takeoff,cutback,climb1 condition)
      V2 = Vstall * 1.2
      M2 = Mstall * 1.2
      CL2 = CLmax / 1.2^2

      ip = iptakeoff
      pare[ieu0, ip] = V2
      pare[ieM0, ip] = M2
      para[iaMach, ip] = M2
      para[iaReunit, ip] = V2 * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2

      #---- set pitch trim by adjusting CLh
      Wf = WTO - Wzero
      rfuel = Wf / parg[igWfuel]
      opt_trim_var = TrimVar.CLHtail
      balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

      CLh2 = para[iaCLh, ip]
      xCG2 = para[iaxCG, ip]
      xCP2 = para[iaxCP, ip]
      xNP2 = para[iaxNP, ip]

      para[iaxCG, 1:iprotate] .= xCG2
      para[iaxCP, 1:iprotate] .= xCP2
      para[iaxNP, 1:iprotate] .= xNP2
      para[iaCLh, 1:iprotate] .= 0.0

      #---- initial guesses for climb points
      ip = ipcutback
      pare[ieu0, ip] = V2
      pare[ieM0, ip] = M2
      para[iaMach, ip] = M2
      para[iaReunit, ip] = V2 * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2
      para[iaxCG, ip] = xCG2
      para[iaxCP, ip] = xCP2
      para[iaxNP, ip] = xNP2
      para[iaCLh, ip] = CLh2

      ip = ipclimb1
      pare[ieu0, ip] = V2
      pare[ieM0, ip] = M2
      para[iaMach, ip] = M2
      para[iaReunit, ip] = V2 * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2
      para[iaxCG, ip] = xCG2
      para[iaxCP, ip] = xCP2
      para[iaxNP, ip] = xNP2
      para[iaCLh, ip] = CLh2

      #---- also use V2 speed for end-of-descent condition, with weight correction
      ip = ipdescentn
      Vrat = sqrt(para[iafracW, ip] / para[iafracW, ipclimb1])
      pare[ieu0, ip] = V2 * Vrat
      pare[ieM0, ip] = M2 * Vrat
      para[iaMach, ip] = M2 * Vrat
      para[iaReunit, ip] = V2 * Vrat * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2
      #
      # ============================================================================
      #---- set up climb points at equal altitude intervals, from altb to altc
      #-    (takeoff ground temperature is neglected here -- std atmosphere is used)
      altb = para[iaalt, iptakeoff]
      altc = para[iaalt, ipcruise1]
      @inbounds for ip = ipclimb1+1:ipclimbn
            frac = float(ip - ipclimb1) / float(ipclimbn - ipclimb1)
            para[iaalt, ip] = altb * (1.0 - frac) + altc * frac

            alt_m = para[iaalt, ip]

            atmos_state = atmos(alt_m, ΔTatmos)
            T0 = atmos_state.T
            p0 = atmos_state.p
            rho0 = atmos_state.ρ
            a0 = atmos_state.a
            mu0 = atmos_state.μ
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = rho0
            pare[iemu0, ip] = mu0

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
            para[iaWbuoy, ip] = (rhocab - rho0) * gee * parg[igcabVol]
      end

      #---- set climb Tt4's from fractions
      fT1 = parg[igfTt4CL1]
      fTn = parg[igfTt4CLn]
      Tt4TO = pare[ieTt4, iptakeoff]
      Tt4CR = pare[ieTt4, ipcruise1]
      @inbounds for ip = ipclimb1:ipclimbn
            frac = float(ip - ipclimb1) /
                   float(ipclimbn - ipclimb1)
            Tfrac = fT1 * (1.0 - frac) + fTn * frac
            pare[ieTt4, ip] = Tt4TO * (1.0 - Tfrac) + Tt4CR * Tfrac
      end
      #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #---- initial values for range, time, weight fraction
      para[iaRange, iprotate:ipclimb1] .= 0.0
      para[iatime, iprotate:ipclimb1] .= 0.0
      para[iafracW, iprotate:ipclimb1] .= WTO / WMTO
      para[iaWbuoy, iprotate:ipclimb1] .= 0.0

      # Initialize climb integrands
      FoW = zeros(Float64, iptotal)
      FFC = zeros(Float64, iptotal)
      Vgi = zeros(Float64, iptotal)

      # integrate trajectory over climb
      @inbounds for ip = ipclimb1:ipclimbn
            if (Ldebug)
                  printstyled("Climb angle integration - ip = ", ip - ipclimb1 + 1, "\n"; color=:light_green)
            end

            # velocity calculation from CL, Weight, altitude
            W = para[iafracW, ip] * WMTO
            CL = para[iaCL, ip]
            ρ = pare[ierho0, ip]
            μ = pare[iemu0, ip]
            Vsound = pare[iea0, ip]
            cosg = cos(para[iagamV, ip])
            BW = W + para[iaWbuoy, ip]
            Wpay = parg[igWpay]
            neng = parg[igneng]

            V = sqrt(2.0 * BW * cosg / (ρ * S * CL))
            Mach = V / Vsound

            para[iaMach, ip] = Mach
            para[iaReunit, ip] = V * ρ / μ

            pare[ieu0, ip] = V
            pare[ieM0, ip] = Mach

                  # Set pitch trim by adjusting CLh
                  Wf = W - Wzero
                  rfuel = Wf / parg[igWfuel]
                  opt_trim_var = TrimVar.CLHtail
                  balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

                  if (ip == ipclimb1)
                        computes_wing_direct = false #use explicitly specified wing cdf, cdp
                  else
                        computes_wing_direct = true #use airfoil database
                  end
                  aircraft_drag!(ac, imission, ip, computes_wing_direct)

                  eng.enginecalc!(ac, "off_design", imission, ip, initializes_engine)

            Ftotal = pare[ieFe, ip] * parg[igneng]
            TSFC = pare[ieTSFC, ip]
            DoL = para[iaCD, ip] / para[iaCL, ip]

            # Calculate improved flight angle
            ϕ = Ftotal / BW
            sing = (ϕ - DoL * sqrt(1.0 - ϕ^2 + DoL^2)) / (1.0 + DoL^2)
            gamV = asin(sing)
            cosg = sqrt(1.0 - sing^2)
            # gamV = atan(sing, cosg)

            para[iagamV, ip] = gamV
            para[iaROC, ip] = sing * V * 60 / ft_to_m #ft per min

            # Store integrands for range and weight integration using a predictor-corrector scheme
            FoW[ip] = Ftotal / (BW * cosg) - DoL

            mfuel = pare[iemfuel, ip]
            FFC[ip] = mfuel * gee / (W * V * cosg)

            Vgi[ip] = 1.0 / (V * cosg)

            Mach = para[iaMach, ip]
            CL = para[iaCL, ip]
            CD = para[iaCD, ip]
            gamV = para[iagamV, ip]

            if (ip > ipclimb1)
                  # Corrector step
                  dh = para[iaalt, ip] - para[iaalt, ip-1]
                  dVsq = pare[ieu0, ip]^2 - pare[ieu0, ip-1]^2

                  FoWavg = 0.5 * (FoW[ip] + FoW[ip-1])
                  FFCavg = 0.5 * (FFC[ip] + FFC[ip-1])
                  Vgiavg = 0.5 * (Vgi[ip] + Vgi[ip-1])

                  dR = (dh + 0.5 * dVsq / gee) / FoWavg
                  dt = dR * Vgiavg
                  rW = exp(-dR * FFCavg)    #Ratio of weights Wᵢ₊₁/Wᵢ

                  para[iaRange, ip] = para[iaRange, ip-1] + dR
                  para[iatime, ip] = para[iatime, ip-1] + dt
                  para[iafracW, ip] = para[iafracW, ip-1] * rW

                  ifirst = false

            end
            if (ip < ipclimbn)
                  # Predictor integration step, forward Euler
                  if (para[iagamV, ip+1] ≤ 0.0)
                        # if gamV guess is not passed in, use previous point as the guess
                        para[iagamV, ip+1] = para[iagamV, ip]
                  end

                  W = para[iafracW, ip+1] * WMTO   # Initial weight fractions have been set in _size_aircraft!
                  CL = para[iaCL, ip+1]
                  ρ = pare[ierho0, ip+1]
                  cosg = cos(para[iagamV, ip+1])

                  BW = W + para[iaWbuoy, ip+1]

                  V = sqrt(2 * BW * cosg / (ρ * S * CL))
                  pare[ieu0, ip+1] = V
                  dh = para[iaalt, ip+1] - para[iaalt, ip]
                  dVsq = pare[ieu0, ip+1]^2 - pare[ieu0, ip]^2

                  dR = (dh + 0.5 * dVsq / gee) / FoW[ip]
                  dt = dR * Vgi[ip]
                  rW = exp(-dR * FFC[ip])


                  para[iaRange, ip+1] = para[iaRange, ip] + dR
                  para[iatime, ip+1] = para[iatime, ip] + dt
                  para[iafracW, ip+1] = para[iafracW, ip] * rW

            end

      end # done integrating climb

      # First cruise point is last climb point
      para[iaRange, ipcruise1] = para[iaRange, ipclimbn]
      para[iatime, ipcruise1] = para[iatime, ipclimbn]
      para[iafracW, ipcruise1] = para[iafracW, ipclimbn]
      para[iaWbuoy, ipcruise1] = para[iaWbuoy, ipclimbn]

      pare[iePLH2, ipcruise1] = pare[iePLH2, ipclimbn]
      pare[ieyg, ipcruise1] = pare[ieyg, ipclimbn]

      # Cruise start
      ip = ipcruise1
      # Set pitch trim by adjusting CLh
      Wf = para[iafracW, ip] * WMTO - Wzero
      rfuel = Wf / parg[igWfuel]
      opt_trim_var = TrimVar.CLHtail
      balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

      if calculate_cruise #If start of cruise has to be calculated (e.g., in off-design)
            # println("Calculating cruise point")
            # Calculate only if requested since for design mission start of cruise is the des point and ∴ already calcualted 
            # Calculate drag
            computes_wing_direct = true
            aircraft_drag!(ac, imission, ip, computes_wing_direct)
            DoL = para[iaCD, ip] / para[iaCL, ip]
            W = para[iafracW, ip] * WMTO
            BW = W + para[iaWbuoy, ip]
            F = BW * (DoL + para[iagamV, ip])
            pare[ieFe, ip] = F / parg[igneng] #Store required thrust for engine calcs
            eng.enginecalc!(ac, "off_design", imission, ip, initializes_engine)

      end

      # set cruise-climb climb angle, from fuel burn rate and atmospheric dp/dz
      TSFC = pare[ieTSFC, ip]
      V = pare[ieu0, ip]
      p0 = pare[iep0, ip]
      ρ0 = pare[ierho0, ip]
      DoL = para[iaCD, ip] / para[iaCL, ip]
      W = para[iafracW, ip] * WMTO
      BW = W + para[iaWbuoy, ip]

      # Calculate Cruise-climb angle:
      # ------
      # γ ≈ dh/dR = dh/dp × dp/dW × dW/dR
      #           = -1/(ρg) × p/W × -̇mf g/(Vcosγ)
      #        ∴γ ≈ ̇mf*p/(ρWV)
      # For conventional aircraft can represent this as function of TSFC: (see TASOPT docs Eq 438)
      # gamVcr1 = DoL*p0*TSFC/(ρ0*gee*V - p0*TSFC)
      gamVcr1 = DoL * p0 * TSFC / (ρ0 * gee * V - p0 * TSFC)

      para[iagamV, ip] = gamVcr1

      Ftotal = pare[ieFe, ip] * parg[igneng]

      cosg = cos(gamVcr1)

      FoW[ip] = Ftotal / (BW * cosg) - DoL
      mfuel = pare[iemfuel, ip]
      FFC[ip] = mfuel * gee / (W * V * cosg)
      Vgi[ip] = 1.0 / (V * cosg)

      #---- set end-of-cruise point "d" using cruise and descent angles, 
      #-     and remaining range legs
      gamVde1 = parm[imgamVDE1]
      gamVden = parm[imgamVDEn]
      gamVdeb = 0.5 * (gamVde1 + gamVden)
      alte = para[iaalt, ipdescentn]
      dRclimb = para[iaRange, ipclimbn] - para[iaRange, ipclimb1]
      dRcruise = (alte - altc - gamVdeb * (Rangetot - dRclimb)) / (gamVcr1 - gamVdeb)

      altd = altc + gamVcr1 * dRcruise

      # ---- build cruise/deviation profile between ipcruise1 and ipcruise2
      Δh_dev = parm[imDeviationDH]
      R_dev_req = parm[imDeviationDL] # off-altitude cruise range between deviation waypoints 2 and 3
      S_dev_req = parm[imDeviationStart]

      dR_available = max(dRcruise, 0.0)
      eps_range = 1.0e-6
      eps_height = 1.0e-3
      # Mission 1 sizes the aircraft and propulsion system. Deviations are
      # operational maneuvers and are evaluated only in off-design missions.
      deviation_requested =
            imission != 1 && abs(Δh_dev) > eps_height && abs(R_dev_req) > eps_range
      Tt4CR = pare[ieTt4, ipcruise1]
      Tt4_avoidance_limit = pare[ieTt4, ipclimbn]
      Tt4_limit_tol = 1.0e-5 * max(abs(Tt4_avoidance_limit), 1.0)

      function deviation_error(reason)
            throw(
                  ArgumentError(
                        "Deviation infeasible for mission $(imission): $(reason). " *
                        "Requested values: dh = $(Δh_dev / ft_to_m) ft, " *
                        "dl = $(R_dev_req / nmi_to_m) nmi, " *
                        "start = $(S_dev_req / nmi_to_m) nmi.",
                  ),
            )
      end

      function estimate_baseline_state(s_target)
            s_eff = max(s_target, 0.0)
            fracW_est = para[iafracW, ipcruise1] * exp(-s_eff * FFC[ipcruise1])
            time_est = para[iatime, ipcruise1] + s_eff * Vgi[ipcruise1]
            alt_est = altc + gamVcr1 * s_eff
            return fracW_est, time_est, alt_est
      end

      function solve_fixed_gamma_point!(
            ipx,
            γ;
            mach = para[iaMach, ipx],
            update_cl = false,
            computes_wing_direct = true,
            t4_initial = nothing,
            max_coupling_iter = 8,
            initialize_engine = initializes_engine,
            cd_tol = 1.0e-7,
            m2_tol = 1.0e-4,
            fe_tol = 1.0e-5,
      )
            para[iagamV, ipx] = γ

            alt_m = para[iaalt, ipx]
            atmos_state = atmos(alt_m, ΔTatmos)
            T0 = atmos_state.T
            p0 = atmos_state.p
            ρ0 = atmos_state.ρ
            a0 = atmos_state.a
            μ0 = atmos_state.μ
            pare[iep0, ipx] = p0
            pare[ieT0, ipx] = T0
            pare[iea0, ipx] = a0
            pare[ierho0, ipx] = ρ0
            pare[iemu0, ipx] = μ0

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
            para[iaWbuoy, ipx] = (rhocab - ρ0) * gee * parg[igcabVol]

            W = para[iafracW, ipx] * WMTO
            BW = W + para[iaWbuoy, ipx]
            cosg = cos(γ)

            V = mach * a0
            para[iaMach, ipx] = mach
            pare[ieM0, ipx] = mach
            pare[ieu0, ipx] = V
            para[iaReunit, ipx] = V * ρ0 / μ0
            if update_cl
                  para[iaCL, ipx] = 2.0 * BW * cosg / (ρ0 * V^2 * S)
            end

            Wf = W - Wzero
            rfuel = Wf / parg[igWfuel]
            opt_trim_var = TrimVar.CLHtail
            balance_aircraft!(ac, imission, ipx, rfuel, rpay, ξpay, opt_trim_var;
                              Ldebug = Ldebug)

            last_cd = Inf
            last_fe = Inf
            for iter = 1:max_coupling_iter
                  m2_in = pare[ieM2, ipx]
                  if iter==1 || !(eng.model isa engine.TurbofanModel)
                        aircraft_drag!(ac, imission, ipx, computes_wing_direct)
                  else
                        aerodynamics.refresh_nacelle_drag!(ac,imission,ipx)
                  end

                  DoL = para[iaCD, ipx] / para[iaCL, ipx]
                  F_required = BW * (DoL * cosg + sin(γ))
                  pare[ieFe, ipx] = F_required / parg[igneng]
                  if iter == 1 && t4_initial !== nothing && eng.model isa engine.TurbofanModel
                        pare[ieTt4, ipx] = t4_initial
                  end

                  if initialize_engine && iter == 1
                        # The caller has copied the converged cruise engine state.
                        # Try that state before discarding it for design-point guesses.
                        # Failed warm solves retain the original cold fallback.
                        engine_initial = copy(view(pare,:,ipx))
                        warm_ok = false
                        try
                              eng.enginecalc!(ac,"off_design",imission,ipx,false)
                              warm_ok = pare[ieConvFail,ipx] == 0
                        catch err
                              err isa Union{InterruptException,MethodError,UndefVarError,BoundsError,TypeError} && rethrow()
                        end
                        if !warm_ok
                              copyto!(view(pare,:,ipx),engine_initial)
                              eng.enginecalc!(ac,"off_design",imission,ipx,true)
                        end
                  else
                        eng.enginecalc!(ac,"off_design",imission,ipx,false)
                  end
                  pare[ieConvFail,ipx] == 0 || deviation_error("engine failed at local state")

                  pare[ieConvFail,ipx] == 0 || deviation_error("local engine solve did not converge")
                  cd = para[iaCD, ipx]
                  fe = pare[ieFe, ipx]
                  m2_out = pare[ieM2, ipx]
                  if (iter == 1 && abs(m2_out-m2_in) <= 0.1*m2_tol) ||
                     (abs(m2_out - m2_in) <= m2_tol &&
                      abs(cd - last_cd) <= cd_tol &&
                      abs(fe - last_fe) <= fe_tol * max(abs(fe), 1.0))
                        break
                  end
                  iter == max_coupling_iter && deviation_error("local drag/engine coupling did not converge")
                  last_cd = cd
                  last_fe = fe
            end

            DoL = para[iaCD, ipx] / para[iaCL, ipx]
            Ftotal = pare[ieFe, ipx] * parg[igneng]
            FoW[ipx] = Ftotal / (BW * cosg) - DoL
            mfuel = pare[iemfuel, ipx]
            FFC[ipx] = mfuel * gee / (W * V * cosg)
            Vgi[ipx] = 1.0 / (V * cosg)
            para[iaROC, ipx] = sin(γ) * V * 60.0 / ft_to_m

            return nothing
      end

      function solve_fixed_Tt4_point!(
            ipx,
            γ,
            Tt4_set;
            mach = para[iaMach, ipx],
            max_coupling_iter = 8,
            initialize_engine = initializes_engine,
            cd_tol = 1.0e-7,
            m2_tol = 1.0e-4,
            fe_tol = 1.0e-5,
      )
            para[iagamV, ipx] = γ

            alt_m = para[iaalt, ipx]
            atmos_state = atmos(alt_m, ΔTatmos)
            T0 = atmos_state.T
            p0 = atmos_state.p
            ρ0 = atmos_state.ρ
            a0 = atmos_state.a
            μ0 = atmos_state.μ
            pare[iep0, ipx] = p0
            pare[ieT0, ipx] = T0
            pare[iea0, ipx] = a0
            pare[ierho0, ipx] = ρ0
            pare[iemu0, ipx] = μ0

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
            para[iaWbuoy, ipx] = (rhocab - ρ0) * gee * parg[igcabVol]

            W = para[iafracW, ipx] * WMTO
            BW = W + para[iaWbuoy, ipx]
            cosg = cos(γ)
            V = mach * a0
            para[iaMach, ipx] = mach
            pare[ieM0, ipx] = mach
            pare[ieu0, ipx] = V
            para[iaReunit, ipx] = V * ρ0 / μ0
            para[iaCL, ipx] = 2.0 * BW * cosg / (ρ0 * V^2 * S)

            Wf = W - Wzero
            rfuel = Wf / parg[igWfuel]
            balance_aircraft!(
                  ac,
                  imission,
                  ipx,
                  rfuel,
                  rpay,
                  ξpay,
                  TrimVar.CLHtail;
                  Ldebug = Ldebug,
            )

            last_cd = Inf
            last_fe = Inf
            for iter = 1:max_coupling_iter
                  m2_in = pare[ieM2, ipx]
                  if iter==1
                        aircraft_drag!(ac, imission, ipx, true)
                  else
                        aerodynamics.refresh_nacelle_drag!(ac,imission,ipx)
                  end
                  pare[ieTt4, ipx] = Tt4_set
                  engine.tfcalc!(
                        wing,
                        eng,
                        parg,
                        view(para, :, ipx),
                        view(pare, :, ipx),
                        ipx,
                        options.ifuel,
                        engine.CalcMode.FixedTt4OffDes,
                        engine.CoolingOpt.FixedCoolingFlowRatio,
                        initialize_engine && iter == 1,
                  )

                  pare[ieConvFail,ipx] == 0 || deviation_error("local engine solve did not converge")
                  cd = para[iaCD, ipx]
                  fe = pare[ieFe, ipx]
                  m2_out = pare[ieM2, ipx]
                  if (iter == 1 && abs(m2_out-m2_in) <= 0.1*m2_tol) ||
                     (abs(m2_out - m2_in) <= m2_tol &&
                      abs(cd - last_cd) <= cd_tol &&
                      abs(fe - last_fe) <= fe_tol * max(abs(fe), 1.0))
                        break
                  end
                  iter == max_coupling_iter && deviation_error("local drag/engine coupling did not converge")
                  last_cd = cd
                  last_fe = fe
            end

            Ftotal = pare[ieFe, ipx] * parg[igneng]
            DoL = para[iaCD, ipx] / para[iaCL, ipx]
            FoW[ipx] = Ftotal / (BW * cosg) - DoL
            mfuel = pare[iemfuel, ipx]
            FFC[ipx] = mfuel * gee / (W * V * cosg)
            Vgi[ipx] = 1.0 / (V * cosg)
            para[iaROC, ipx] = sin(γ) * V * 60.0 / ft_to_m
            return Ftotal
      end

      function climb_gamma_from_fixed_Tt4_endpoint!(ipx, γ, Tt4_set)
            Ftotal = solve_fixed_Tt4_point!(
                  ipx,
                  γ,
                  Tt4_set;
                  initialize_engine = false,
            )
            W = para[iafracW, ipx] * WMTO
            BW = W + para[iaWbuoy, ipx]
            DoL = para[iaCD, ipx] / para[iaCL, ipx]
            ϕ = Ftotal / BW
            radicand = 1.0 - ϕ^2 + DoL^2
            radicand >= 0.0 || return NaN
            sing = (ϕ - DoL * sqrt(radicand)) / (1.0 + DoL^2)
            return asin(clamp(sing, -1.0, 1.0))
      end

      segment_cache=Dict{NTuple{9,Float64},MissionSegmentCacheEntry}()
      function propagate_segment!(ipx, s0, s1, h0, h1, w0, t0, gamma, mach;
                                  enforce_limit=false)
            key=(Float64(ipx),s0,s1,h0,h1,w0,t0,gamma,mach)
            if integration_cache && haskey(segment_cache,key)
                  entry=segment_cache[key]
                  if enforce_limit && !allow_deviation_t4_violation &&
                     entry.maximum_temperature > Tt4_avoidance_limit+Tt4_limit_tol
                        deviation_error("cached converged segment exceeds climb temperature limit")
                  end
                  @views para[:,ipx].=entry.aero
                  @views pare[:,ipx].=entry.engine
                  FFC[ipx],Vgi[ipx],FoW[ipx]=entry.rates
                  return entry.result
            end
            temperatures = fill(NaN,4)
            accepted_max_temperature=-Inf
            @views para[:,ipx] .= para[:,ipcruise1]
            @views pare[:,ipx] .= pare[:,ipcruise1]
            first_state=true
            function rates!(s,w)
                  fraction = s1 == s0 ? 0.0 : (s-s0)/(s1-s0)
                  para[iaRange,ipx] = para[iaRange,ipcruise1] + s
                  para[iaalt,ipx] = h0 + fraction*(h1-h0)
                  para[iafracW,ipx] = w
                  solve_fixed_gamma_point!(ipx,gamma;mach,update_cl=true,
                      t4_initial=first_state ? Tt4CR : nothing,max_coupling_iter=8,initialize_engine=first_state)
                  first_state && (accepted_max_temperature=pare[ieTt4,ipx])
                  first_state=false
                  temperatures[1]=temperatures[2];temperatures[2]=temperatures[3];temperatures[3]=temperatures[4];temperatures[4]=pare[ieTt4,ipx]
                  return FFC[ipx],Vgi[ipx]
            end
            function validate!()
                  accepted_max_temperature=max(accepted_max_temperature,maximum(temperatures))
                  if enforce_limit && !allow_deviation_t4_violation &&
                     accepted_max_temperature > Tt4_avoidance_limit+Tt4_limit_tol
                        deviation_error("converged segment state exceeds climb temperature limit")
                  end
            end
            result=integrate_mission_segment!(rates!,s0,s1,w0,t0;
                       rtol=integration_rtol,max_step=min(integration_max_step, (s1>s0 && abs(h1-h0)>1e-9) ? (s1-s0)*1000*ft_to_m/abs(h1-h0) : integration_max_step),validate!)
            para[iatime,ipx]=result.time
            if integration_cache
                  segment_cache[key]=MissionSegmentCacheEntry(result,copy(view(para,:,ipx)),
                      copy(view(pare,:,ipx)),(FFC[ipx],Vgi[ipx],FoW[ipx]),accepted_max_temperature)
            end
            return result
      end

      use_deviation = false
      if deviation_requested
            if abs(Δh_dev) <= eps_height
                  deviation_error("a nonzero deviation length requires a nonzero deviation height")
            elseif R_dev_req <= eps_range
                  deviation_error("a nonzero deviation height requires a positive deviation length")
            elseif S_dev_req < -eps_range
                  deviation_error("deviation start cannot be negative")
            elseif dR_available <= eps_range
                  deviation_error("no cruise span is available for a deviation")
            end

            if R_dev_req > dR_available + eps_range
                  deviation_error(
                        "requested deviation cruise length exceeds available cruise span of $(dR_available / nmi_to_m) nmi",
                  )
            elseif S_dev_req + R_dev_req > dR_available + eps_range
                  deviation_error(
                        "requested deviation start plus deviation cruise length exceeds the available cruise span before accounting for transition legs",
                  )
            end

            # The normal cruise Mach is prescribed at the ordinary cruise
            # endpoint.  A distinct avoidance Mach can be prescribed in the
            # first deviation slot before the mission solve.  Existing models
            # remain backward compatible because those slots already contain
            # the normal cruise Mach.
            Mach_normal = para[iaMach, ipcruise1]
            Mach_avoidance_requested = para[iaMach, ipdeviation1]
            Mach_avoidance =
                  isfinite(Mach_avoidance_requested) && Mach_avoidance_requested > 0.0 ?
                  Mach_avoidance_requested : Mach_normal
            range_start = para[iaRange, ipcruise1]

            D1 = 1
            D2 = 2
            D3 = 3
            D4 = 4
            D5 = 5
            D6 = 6
            D7 = 7
            D8 = 8
            ndev = 8

            dev_s = zeros(Float64, ndev)
            dev_offset = zeros(Float64, ndev)
            dev_alt = zeros(Float64, ndev)
            dev_time = zeros(Float64, ndev)
            dev_fracW = zeros(Float64, ndev)
            dev_gamma = fill(gamVcr1, ndev)
            dev_mode = fill(:cruise, ndev)
            # D1 and D8 are the normal-cruise states immediately before and
            # after the maneuver. D2:D7 span the avoidance climb, hold, and
            # descent, including the colocated phase-switch states.
            dev_mach = Float64[
                  Mach_normal,
                  Mach_avoidance,
                  Mach_avoidance,
                  Mach_avoidance,
                  Mach_avoidance,
                  Mach_avoidance,
                  Mach_avoidance,
                  Mach_normal,
            ]
            dev_para = Matrix{Float64}(undef, size(para, 1), ndev)
            dev_pare = Matrix{Float64}(undef, size(pare, 1), ndev)
            dev_FFC = zeros(Float64, ndev)
            dev_Vgi = zeros(Float64, ndev)
            dev_FoW = zeros(Float64, ndev)

            return_gamma_mode = lowercase(strip(get(ENV, "TASOPT_AVOIDANCE_RETURN_GAMMA_MODE", "current")))
            if !(return_gamma_mode in ("current", "default", "mission", "symmetric", "mirror"))
                  deviation_error(
                        "unknown TASOPT_AVOIDANCE_RETURN_GAMMA_MODE=$(return_gamma_mode); use current or symmetric",
                  )
            end
            symmetric_return_gamma = return_gamma_mode in ("symmetric", "mirror")

            function avoidance_descent_gamma(γ_climb)
                  return symmetric_return_gamma ? -abs(γ_climb) : parm[imgamVDE1]
            end

            function build_eight_point_deviation_profile(γ_climb, γ_descent)
                  γ_first = Δh_dev > 0.0 ? γ_climb : γ_descent
                  γ_return = Δh_dev > 0.0 ? γ_descent : γ_climb
                  den_first = γ_first - gamVcr1
                  den_return = γ_return - gamVcr1
                  eps_den = 1.0e-8

                  if abs(den_first) < eps_den
                        deviation_error("first deviation transition angle matches the baseline cruise angle, so D3/D4 cannot leave the baseline path")
                  elseif abs(den_return) < eps_den
                        deviation_error("return deviation transition angle matches the baseline cruise angle, so D7/D8 cannot return to the baseline path")
                  end

                  Δs_first = Δh_dev / γ_first
                  if Δs_first <= eps_range
                        deviation_error("computed first transition leg length is not positive")
                  end

                  offset_hold = (γ_first - gamVcr1) * Δs_first
                  if abs(offset_hold) < eps_height
                        deviation_error("the requested deviation height stays on the baseline path within tolerance")
                  end

                  Δs_return = -offset_hold / den_return
                  if Δs_return <= eps_range
                        deviation_error("computed return transition leg length is not positive")
                  end

                  total_length = Δs_first + R_dev_req + Δs_return
                  if total_length > dR_available + eps_range
                        deviation_error(
                              "minimum achievable deviation span is $(total_length / nmi_to_m) nmi, which exceeds the available cruise span of $(dR_available / nmi_to_m) nmi",
                        )
                  elseif S_dev_req + total_length > dR_available + eps_range
                        deviation_error(
                              "requested deviation start plus transition legs and deviation cruise length exceeds the available cruise span",
                        )
                  end

                  return (
                        γ_first = γ_first,
                        γ_return = γ_return,
                        Δs_first = Δs_first,
                        Δs_return = Δs_return,
                        R_hold = R_dev_req,
                        total_length = total_length,
                  )
            end

            function apply_eight_point_geometry!(profile, γ_climb, γ_descent)
                  s12 = S_dev_req
                  s34 = s12 + profile.Δs_first
                  s56 = s34 + profile.R_hold
                  s78 = s56 + profile.Δs_return

                  dev_s[D1] = s12
                  dev_s[D2] = s12
                  dev_s[D3] = s34
                  dev_s[D4] = s34
                  dev_s[D5] = s56
                  dev_s[D6] = s56
                  dev_s[D7] = s78
                  dev_s[D8] = s78

                  dev_offset[D1] = 0.0
                  dev_offset[D2] = 0.0
                  dev_offset[D3] = profile.γ_first * profile.Δs_first - gamVcr1 * profile.Δs_first
                  dev_offset[D4] = dev_offset[D3]
                  dev_offset[D5] = dev_offset[D3]
                  dev_offset[D6] = dev_offset[D3]
                  dev_offset[D7] = 0.0
                  dev_offset[D8] = 0.0

                  @inbounds for i = 1:ndev
                        dev_alt[i] = altc + gamVcr1 * dev_s[i] + dev_offset[i]
                        fracW_est, time_est, _ = estimate_baseline_state(dev_s[i])
                        dev_fracW[i] = fracW_est
                        dev_time[i] = time_est
                  end

                  if Δh_dev > 0.0
                        dev_gamma[D1] = gamVcr1
                        dev_gamma[D2] = γ_climb
                        dev_gamma[D3] = γ_climb
                        dev_gamma[D4] = gamVcr1
                        dev_gamma[D5] = gamVcr1
                        dev_gamma[D6] = γ_descent
                        dev_gamma[D7] = γ_descent
                        dev_gamma[D8] = gamVcr1

                        dev_mode[D1] = :cruise
                        dev_mode[D2] = :climb
                        dev_mode[D3] = :climb
                        dev_mode[D4] = :cruise
                        dev_mode[D5] = :cruise
                        dev_mode[D6] = :descent
                        dev_mode[D7] = :descent
                        dev_mode[D8] = :cruise
                  else
                        dev_gamma[D1] = gamVcr1
                        dev_gamma[D2] = γ_descent
                        dev_gamma[D3] = γ_descent
                        dev_gamma[D4] = gamVcr1
                        dev_gamma[D5] = gamVcr1
                        dev_gamma[D6] = γ_climb
                        dev_gamma[D7] = γ_climb
                        dev_gamma[D8] = gamVcr1

                        dev_mode[D1] = :cruise
                        dev_mode[D2] = :descent
                        dev_mode[D3] = :descent
                        dev_mode[D4] = :cruise
                        dev_mode[D5] = :cruise
                        dev_mode[D6] = :climb
                        dev_mode[D7] = :climb
                        dev_mode[D8] = :cruise
                  end
            end

            function populate_deviation_point!(
                  ipx,
                  γ,
                  mode;
                  mach = Mach_normal,
                  enforce_climb_Tt4_limit = true,
                  reason_label = "$(mode) deviation point",
                  couple_engine = false,
            )
                  couple_for_limit =
                        couple_engine ||
                        (
                              enforce_climb_Tt4_limit &&
                              mode == :climb &&
                              eng.model isa engine.TurbofanModel
                        )

                  try
                        solve_fixed_gamma_point!(
                              ipx,
                              γ;
                              mach,
                              update_cl = true,
                              computes_wing_direct = true,
                              t4_initial = Tt4CR,
                              max_coupling_iter = 8,
                        )
                  catch err
                        err isa Union{InterruptException,MethodError,UndefVarError,BoundsError,TypeError} && rethrow()
                        deviation_error("required-thrust engine solve failed for the $(reason_label): $(sprint(showerror,err))")
                  end

                  if enforce_climb_Tt4_limit &&
                     mode == :climb &&
                     eng.model isa engine.TurbofanModel &&
                     !allow_deviation_t4_violation &&
                     pare[ieTt4, ipx] > Tt4_avoidance_limit + Tt4_limit_tol
                        deviation_error(
                              "required Tt4 $(pare[ieTt4, ipx]) K exceeds avoidance/climb Tt4 limit $(Tt4_avoidance_limit) K for the $(reason_label)",
                        )
                  end
            end

            function solve_local_deviation_state!(
                  idx;
                  enforce_climb_Tt4_limit = true,
                  couple_engine = false,
            )
                  @views para[:, iptest] .= para[:, ipcruise1]
                  @views pare[:, iptest] .= pare[:, ipcruise1]
                  para[iaRange, iptest] = range_start + dev_s[idx]
                  para[iaalt, iptest] = dev_alt[idx]
                  para[iatime, iptest] = dev_time[idx]
                  para[iafracW, iptest] = dev_fracW[idx]
                  para[iaCL, iptest] = para[iaCL, ipcruise1]

                  populate_deviation_point!(
                        iptest,
                        dev_gamma[idx],
                        dev_mode[idx];
                        mach = dev_mach[idx],
                        enforce_climb_Tt4_limit = enforce_climb_Tt4_limit,
                        reason_label = "D$(idx) $(dev_mode[idx]) deviation point",
                        couple_engine = couple_engine,
                  )

                  @views dev_para[:, idx] .= para[:, iptest]
                  @views dev_pare[:, idx] .= pare[:, iptest]
                  dev_FoW[idx] = FoW[iptest]
                  dev_FFC[idx] = FFC[iptest]
                  dev_Vgi[idx] = Vgi[iptest]
            end

            function seed_local_from_global!(ip_start, idx_end)
                  s0=para[iaRange,ip_start]-range_start
                  result=propagate_segment!(iptest,s0,dev_s[idx_end],para[iaalt,ip_start],dev_alt[idx_end],
                      para[iafracW,ip_start],para[iatime,ip_start],gamVcr1,Mach_normal)
                  dev_time[idx_end]=result.time;dev_fracW[idx_end]=result.weight
                  @views dev_para[:,idx_end].=para[:,iptest]
                  @views dev_pare[:,idx_end].=pare[:,iptest]
                  dev_FoW[idx_end]=FoW[iptest];dev_FFC[idx_end]=FFC[iptest];dev_Vgi[idx_end]=Vgi[iptest]
            end

            function integrate_local_to_local!(idx_start,idx_end;
                    enforce_climb_Tt4_limit=true,couple_engine=false)
                  result=propagate_segment!(iptest,dev_s[idx_start],dev_s[idx_end],
                      dev_alt[idx_start],dev_alt[idx_end],dev_fracW[idx_start],dev_time[idx_start],
                      dev_gamma[idx_end],dev_mach[idx_end];
                      enforce_limit=enforce_climb_Tt4_limit && dev_mode[idx_end]==:climb)
                  dev_time[idx_end]=result.time;dev_fracW[idx_end]=result.weight
                  @views dev_para[:,idx_end].=para[:,iptest]
                  @views dev_pare[:,idx_end].=pare[:,iptest]
                  dev_FoW[idx_end]=FoW[iptest];dev_FFC[idx_end]=FFC[iptest];dev_Vgi[idx_end]=Vgi[iptest]
            end

            function switch_local_phase!(
                  idx_start,
                  idx_end;
                  enforce_climb_Tt4_limit = true,
                  couple_engine = false,
            )
                  dev_time[idx_end] = dev_time[idx_start]
                  dev_fracW[idx_end] = dev_fracW[idx_start]
                  solve_local_deviation_state!(
                        idx_end;
                        enforce_climb_Tt4_limit = enforce_climb_Tt4_limit,
                        couple_engine = couple_engine,
                  )
            end

            function solve_eight_point_deviation_states!(
                  profile,
                  γ_climb,
                  γ_descent;
                  enforce_climb_Tt4_limit = true,
                  couple_engine = false,
            )
                  apply_eight_point_geometry!(profile, γ_climb, γ_descent)

                  seed_local_from_global!(ipcruise1, D1)
                  switch_local_phase!(
                        D1,
                        D2;
                        enforce_climb_Tt4_limit = enforce_climb_Tt4_limit,
                        couple_engine = couple_engine,
                  )
                  integrate_local_to_local!(
                        D2,
                        D3;
                        enforce_climb_Tt4_limit = enforce_climb_Tt4_limit,
                        couple_engine = couple_engine,
                  )
                  switch_local_phase!(
                        D3,
                        D4;
                        enforce_climb_Tt4_limit = false,
                        couple_engine = couple_engine,
                  )
                  integrate_local_to_local!(
                        D4,
                        D5;
                        enforce_climb_Tt4_limit = false,
                        couple_engine = couple_engine,
                  )
                  switch_local_phase!(
                        D5,
                        D6;
                        enforce_climb_Tt4_limit = enforce_climb_Tt4_limit,
                        couple_engine = couple_engine,
                  )
                  integrate_local_to_local!(
                        D6,
                        D7;
                        enforce_climb_Tt4_limit = enforce_climb_Tt4_limit,
                        couple_engine = couple_engine,
                  )
                  switch_local_phase!(
                        D7,
                        D8;
                        enforce_climb_Tt4_limit = false,
                        couple_engine = couple_engine,
                  )
            end

            function solve_active_climb_margin_states!(
                  profile,
                  γ_climb,
                  γ_descent;
                  couple_engine = false,
            )
                  apply_eight_point_geometry!(profile, γ_climb, γ_descent)

                  seed_local_from_global!(ipcruise1, D1)

                  if Δh_dev > 0.0
                        switch_local_phase!(
                              D1,
                              D2;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                        integrate_local_to_local!(
                              D2,
                              D3;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                  else
                        switch_local_phase!(
                              D1,
                              D2;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                        integrate_local_to_local!(
                              D2,
                              D3;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                        switch_local_phase!(
                              D3,
                              D4;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                        integrate_local_to_local!(
                              D4,
                              D5;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                        switch_local_phase!(
                              D5,
                              D6;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                        integrate_local_to_local!(
                              D6,
                              D7;
                              enforce_climb_Tt4_limit = false,
                              couple_engine = couple_engine,
                        )
                  end
            end

            function active_climb_state_indices()
                  return Δh_dev > 0.0 ? (D2, D3) : (D6, D7)
            end

            function active_climb_endpoint_gamma_capability(γ_climb)
                  γ_descent = avoidance_descent_gamma(γ_climb)
                  profile = build_eight_point_deviation_profile(γ_climb, γ_descent)
                  solve_active_climb_margin_states!(
                        profile,
                        γ_climb,
                        γ_descent;
                        couple_engine = true,
                  )

                  idx_start, idx_end = active_climb_state_indices()
                  γ_capability = Inf
                  for idx in (idx_start, idx_end)
                        @views para[:, iptest] .= dev_para[:, idx]
                        @views pare[:, iptest] .= dev_pare[:, idx]
                        γ_endpoint = climb_gamma_from_fixed_Tt4_endpoint!(
                              iptest,
                              γ_climb,
                              Tt4_avoidance_limit,
                        )
                        isfinite(γ_endpoint) || return -Inf
                        γ_capability = min(γ_capability, γ_endpoint)
                  end
                  return γ_capability
            end

            function active_climb_endpoint_Tt4_required(γ_climb)
                  γ_descent = avoidance_descent_gamma(γ_climb)
                  profile = build_eight_point_deviation_profile(γ_climb, γ_descent)
                  solve_active_climb_margin_states!(
                        profile,
                        γ_climb,
                        γ_descent;
                        couple_engine = true,
                  )
                  idx_start, idx_end = active_climb_state_indices()
                  return max(
                        dev_pare[ieTt4, idx_start],
                        dev_pare[ieTt4, idx_end],
                  )
            end

            function minimum_climb_gamma_for_geometry()
                  transition_range_available =
                        dR_available - S_dev_req - R_dev_req
                  transition_range_available > eps_range || deviation_error(
                        "no cruise range remains for the deviation transition legs",
                  )

                  γ_floor = gamVcr1 + 1.0e-7
                  γ_geom = γ_floor
                  if symmetric_return_gamma
                        height = abs(Δh_dev)
                        γ_geom =
                              Δh_dev > 0.0 ?
                              2.0 * height / transition_range_available - gamVcr1 :
                              2.0 * height / transition_range_available + gamVcr1
                  else
                        γ_descent = parm[imgamVDE1]
                        γ_descent < gamVcr1 - 1.0e-8 || deviation_error(
                              "the prescribed deviation descent angle must be below the baseline cruise angle",
                        )

                        if Δh_dev > 0.0
                              range_ratio =
                                    transition_range_available *
                                    (gamVcr1 - γ_descent) / Δh_dev
                              range_ratio > 1.0 || deviation_error(
                                    "minimum achievable deviation span exceeds the available cruise span",
                              )
                              γ_geom = γ_descent / (1.0 - range_ratio)
                        else
                              first_leg_range = Δh_dev / γ_descent
                              remaining_return_range =
                                    transition_range_available - first_leg_range
                              remaining_return_range > eps_range || deviation_error(
                                    "minimum achievable deviation span exceeds the available cruise span",
                              )
                              offset_after_first_leg =
                                    (γ_descent - gamVcr1) * first_leg_range
                              γ_geom =
                                    gamVcr1 -
                                    offset_after_first_leg / remaining_return_range
                        end
                  end

                  isfinite(γ_geom) || deviation_error(
                        "could not determine the minimum climb angle required by the deviation geometry",
                  )
                  # Move an infinitesimal distance into the feasible side so the
                  # range check is not decided by roundoff at exact equality.
                  return max(γ_geom, γ_floor) + 1.0e-10
            end

            gamma_geometry_min_diagnostic = NaN
            gamma_thrust_limit_diagnostic = NaN
            t4_geometry_required_diagnostic = NaN

            function choose_climb_gamma()
                  γ_geom_min = minimum_climb_gamma_for_geometry()
                  gamma_geometry_min_diagnostic = γ_geom_min

                  if !(eng.model isa engine.TurbofanModel)
                        γ_nominal_climb = maximum(view(para, iagamV, ipclimb1:ipclimbn))
                        if !isfinite(γ_nominal_climb) || γ_nominal_climb <= gamVcr1
                              γ_nominal_climb = γ_geom_min
                        end
                        return max(γ_nominal_climb, γ_geom_min)
                  end

                  # The geometry-limited angle is a benign off-design state and
                  # gives the fixed-Tt4 solve a reliable initial engine state.
                  γ_climb = γ_geom_min
                  # A cheap baseline-weight prediction supplies only an initial
                  # angle. The accepted angle still uses converged local states
                  # and the unchanged strict fixed-T4/fixed-thrust checks below.
                  try
                        for prediction in 1:3
                              γ_descent = avoidance_descent_gamma(γ_climb)
                              profile = build_eight_point_deviation_profile(γ_climb,γ_descent)
                              apply_eight_point_geometry!(profile,γ_climb,γ_descent)
                              capability=Inf
                              for idx in active_climb_state_indices()
                                    @views para[:,iptest].=para[:,ipcruise1]
                                    @views pare[:,iptest].=pare[:,ipcruise1]
                                    para[iaalt,iptest]=dev_alt[idx]
                                    para[iafracW,iptest]=dev_fracW[idx]
                                    para[iaMach,iptest]=dev_mach[idx]
                                    capability=min(capability,climb_gamma_from_fixed_Tt4_endpoint!(iptest,γ_climb,Tt4_avoidance_limit))
                              end
                              isfinite(capability) || break
                              γ_climb=max(γ_geom_min,capability)
                        end
                  catch err
                        err isa ArgumentError || rethrow()
                        γ_climb=γ_geom_min
                  end

                  # Endpoint weight, altitude, and drag depend weakly on the
                  # prescribed angle. Iterate only that self-consistency: each
                  # pass evaluates both endpoints once at the fixed Tt4 limit.
                  γ_tol = 1.0e-7
                  γ_direct = NaN
                  for _ = 1:12
                        γ_thrust_limit =
                              active_climb_endpoint_gamma_capability(γ_climb)
                        gamma_thrust_limit_diagnostic = γ_thrust_limit
                        if !isfinite(γ_thrust_limit) ||
                           γ_thrust_limit < γ_geom_min - γ_tol
                              if allow_deviation_t4_violation && isfinite(γ_thrust_limit)
                                    t4_geometry_required_diagnostic =
                                          active_climb_endpoint_Tt4_required(γ_geom_min)
                                    return γ_geom_min
                              end
                              deviation_error(
                                    "not enough thrust available under the avoidance/climb Tt4 limit for the deviation climb " *
                                    "(minimum geometry gamma=$(γ_geom_min), " *
                                    "Tt4-limited gamma=$(γ_thrust_limit), " *
                                    "gamma margin=$(γ_geom_min - γ_thrust_limit))",
                              )
                        end

                        if abs(γ_thrust_limit - γ_climb) <= γ_tol
                              γ_direct = max(γ_thrust_limit, γ_geom_min)
                              break
                        end
                        γ_climb = γ_thrust_limit
                  end

                  isfinite(γ_direct) || deviation_error(
                        "the Tt4-limited deviation climb angle did not converge",
                  )

                  # The fixed-Tt4 result is the physical angle prediction. The
                  # mission point is stored by the fixed-thrust off-design
                  # solver, so correct their small numerical inverse mismatch.
                  Tt4_direct = active_climb_endpoint_Tt4_required(γ_direct)
                  if Tt4_direct <= Tt4_avoidance_limit + Tt4_limit_tol
                        return γ_direct
                  end

                  γ_low = γ_geom_min
                  Tt4_low = active_climb_endpoint_Tt4_required(γ_low)
                  t4_geometry_required_diagnostic = Tt4_low
                  if !isfinite(Tt4_low) || Tt4_low > Tt4_avoidance_limit + Tt4_limit_tol
                        if allow_deviation_t4_violation && isfinite(Tt4_low)
                              return γ_low
                        end
                        deviation_error(
                              "not enough thrust available under the avoidance/climb Tt4 limit for the deviation climb " *
                              "(geometry-required Tt4=$(Tt4_low) K, " *
                              "limit=$(Tt4_avoidance_limit) K, " *
                              "limit tolerance=$(Tt4_limit_tol) K, " *
                              "Tt4 margin=$(Tt4_low - (Tt4_avoidance_limit + Tt4_limit_tol)) K)",
                        )
                  end

                  γ_high = γ_direct
                  Tt4_high = Tt4_direct
                  for _ = 1:15
                        temperature_span = Tt4_high - Tt4_low
                        temperature_span > 0.0 || break
                        γ_trial = γ_low +
                                  (Tt4_avoidance_limit - Tt4_low) *
                                  (γ_high - γ_low) / temperature_span
                        γ_trial = clamp(γ_trial, γ_low, γ_high)
                        Tt4_trial = active_climb_endpoint_Tt4_required(γ_trial)
                        isfinite(Tt4_trial) || break

                        if abs(Tt4_trial - Tt4_avoidance_limit) <= Tt4_limit_tol
                              return γ_trial
                        elseif Tt4_trial < Tt4_avoidance_limit
                              γ_low = γ_trial
                              Tt4_low = Tt4_trial
                        else
                              γ_high = γ_trial
                              Tt4_high = Tt4_trial
                        end
                  end

                  deviation_error(
                        "the fixed-Tt4 and fixed-thrust deviation endpoint solves did not converge consistently",
                  )
            end

            function write_local_state_to_ip!(
                  idx,
                  ipx,
            )
                  @views para[:, ipx] .= dev_para[:, idx]
                  @views pare[:, ipx] .= dev_pare[:, idx]
                  FFC[ipx] = dev_FFC[idx]
                  Vgi[ipx] = dev_Vgi[idx]
                  FoW[ipx] = dev_FoW[idx]
            end

            function solve_final_cruise_point!()
                  propagate_segment!(ipcruisen,dev_s[D8],dRcruise,dev_alt[D8],altd,
                      dev_fracW[D8],dev_time[D8],gamVcr1,Mach_normal)
            end

            para_iptest = copy(view(para, :, iptest))
            pare_iptest = copy(view(pare, :, iptest))
            try
                  γ_climb = choose_climb_gamma()
                  γ_descent = avoidance_descent_gamma(γ_climb)
                  if γ_descent >= gamVcr1 - 1.0e-8
                        deviation_error("avoidance descent angle is not below the baseline cruise angle")
                  end

                  final_profile = build_eight_point_deviation_profile(γ_climb, γ_descent)
                  solve_eight_point_deviation_states!(
                        final_profile,
                        γ_climb,
                        γ_descent;
                        enforce_climb_Tt4_limit = true,
                        couple_engine = false,
                  )

                  for idx in D1:D8
                        write_local_state_to_ip!(idx, ipmaneuver1+idx-1)
                  end
                  solve_final_cruise_point!()

                  achieved_Δh = dev_alt[D3] - dev_alt[D1]
                  achieved_R_dev = dev_s[D5] - dev_s[D3]
                  paired_range_ok =
                        isapprox(dev_s[D1], dev_s[D2]; atol = eps_range, rtol = 0.0) &&
                        isapprox(dev_s[D3], dev_s[D4]; atol = eps_range, rtol = 0.0) &&
                        isapprox(dev_s[D5], dev_s[D6]; atol = eps_range, rtol = 0.0) &&
                        isapprox(dev_s[D7], dev_s[D8]; atol = eps_range, rtol = 0.0)
                  paired_time_ok =
                        isapprox(dev_time[D1], dev_time[D2]; atol = 1.0e-8, rtol = 0.0) &&
                        isapprox(dev_time[D3], dev_time[D4]; atol = 1.0e-8, rtol = 0.0) &&
                        isapprox(dev_time[D5], dev_time[D6]; atol = 1.0e-8, rtol = 0.0) &&
                        isapprox(dev_time[D7], dev_time[D8]; atol = 1.0e-8, rtol = 0.0)

                  achieved_Δh_nodes = para[iaalt, ipdeviation2] - para[iaalt, ipdeviation1]
                  achieved_R_nodes = para[iaRange, ipdeviation3] - para[iaRange, ipdeviation2]
                  achieved_hold_offset_nodes =
                        para[iaalt, ipdeviation3] - para[iaalt, ipdeviation2] -
                        gamVcr1 * achieved_R_nodes

                  if !paired_range_ok || !paired_time_ok
                        deviation_error("internal D1-D8 phase pairs do not share range and time")
                  elseif !isapprox(achieved_Δh, Δh_dev; atol = eps_height, rtol = 0.0)
                        deviation_error(
                              "internal geometry failed to match the requested D1/D2 to D3/D4 altitude change",
                        )
                  elseif !isapprox(achieved_R_dev, R_dev_req; atol = eps_range, rtol = 0.0)
                        deviation_error(
                              "internal geometry failed to match the requested D3/D4 to D5/D6 deviation cruise length",
                        )
                  elseif !isapprox(achieved_Δh_nodes, Δh_dev; atol = eps_height, rtol = 0.0)
                        deviation_error("constructed waypoint altitudes do not satisfy the requested deviation altitude change")
                  elseif !isapprox(achieved_R_nodes, R_dev_req; atol = eps_range, rtol = 0.0)
                        deviation_error("constructed waypoint ranges do not satisfy the requested deviation cruise length")
                  elseif !isapprox(achieved_hold_offset_nodes, 0.0; atol = eps_height, rtol = 0.0)
                        deviation_error("constructed waypoint altitudes do not hold the requested deviation offset during off-altitude cruise")
                  end

                  _store_deviation_diagnostic!(
                        ac,
                        imission,
                        range_start,
                        dev_s,
                        dev_alt,
                        dev_time,
                        dev_fracW,
                        dev_gamma,
                        dev_mode,
                        gamma_geometry_min_diagnostic,
                        gamma_thrust_limit_diagnostic,
                        t4_geometry_required_diagnostic,
                        Tt4_avoidance_limit,
                        Tt4_limit_tol,
                        dev_para,
                        dev_pare,
                        dev_FFC,
                        dev_Vgi,
                        dev_FoW,
                  )

                  use_deviation = true
            finally
                  @views para[:, iptest] .= para_iptest
                  @views pare[:, iptest] .= pare_iptest
            end
      end

      if !use_deviation
            # Same state equations and quadrature as pre/hold/post-deviation cruise.
            result=propagate_segment!(ipcruisen,0.,dRcruise,altc,altd,
                para[iafracW,ipcruise1],para[iatime,ipcruise1],gamVcr1,para[iaMach,ipcruise1])
            dtcruise=result.time-para[iatime,ipcruise1]
            rWcruise=result.weight/para[iafracW,ipcruise1]

            # set intermediate points over cruise, if any, just by interpolating
            for ip = ipcruise1+1:ipcruisen-1
                  frac = float(ip - ipcruise1) / float(ipcruisen - ipcruise1)
                  Mach = para[iaMach, ip]
                  para[iaalt, ip] = altc * (1.0 - frac) + altd * frac
                  alt_m = para[iaalt, ip]
                  atmos_state = atmos(alt_m, ΔTatmos)
                  T0 = atmos_state.T
                  p0 = atmos_state.p
                  rho0 = atmos_state.ρ
                  a0 = atmos_state.a
                  mu0 = atmos_state.μ
                  pare[iep0, ip] = p0
                  pare[ieT0, ip] = T0
                  pare[iea0, ip] = a0
                  pare[ierho0, ip] = rho0
                  pare[iemu0, ip] = mu0
                  pare[ieM0, ip] = Mach
                  pare[ieu0, ip] = Mach * a0
                  para[iaReunit, ip] = Mach * a0 * rho0 / mu0

                  rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
                  para[iaWbuoy, ip] = (rhocab - rho0) * gee * parg[igcabVol]

                  para[iaRange, ip] = para[iaRange, ipcruise1] + dRcruise * frac
                  para[iatime, ip] = para[iatime, ipcruise1] + dtcruise * frac
                  para[iafracW, ip] = para[iafracW, ipcruise1] * rWcruise^frac

            end
      end

      # Descent
      ip = ipdescent1
      pare[iep0, ip] = pare[iep0, ipcruisen]
      pare[ieT0, ip] = pare[ieT0, ipcruisen]
      pare[iea0, ip] = pare[iea0, ipcruisen]
      pare[ierho0, ip] = pare[ierho0, ipcruisen]
      pare[iemu0, ip] = pare[iemu0, ipcruisen]
      pare[ieM0, ip] = pare[ieM0, ipcruisen]
      pare[ieu0, ip] = pare[ieu0, ipcruisen]

      para[iaMach, ip] = para[iaMach, ipcruisen]
      para[iaReunit, ip] = para[iaReunit, ipcruisen]
      para[iaalt, ip] = para[iaalt, ipcruisen]

      para[iaRange, ip] = para[iaRange, ipcruisen]
      para[iatime, ip] = para[iatime, ipcruisen]
      para[iafracW, ip] = para[iafracW, ipcruisen]
      para[iaWbuoy, ip] = para[iaWbuoy, ipcruisen]

      Rd = para[iaRange, ipdescent1]
      Re = para[iaRange, ipdescentn]
      altd = para[iaalt, ipdescent1]
      alte = para[iaalt, ipdescentn]
      para[iagamV, ipdescent1] = gamVde1
      para[iagamV, ipdescentn] = gamVden

      for ip = ipdescent1+1:ipdescentn-1
            frac = float(ip - ipdescent1) / float(ipdescentn - ipdescent1)
            R = Rd * (1.0 - frac) + Re * frac

            alt = altd + (Re - Rd) * (gamVde1 * (frac - 0.5 * frac^2) + gamVden * 0.5 * frac^2)
            gamVde = gamVde1 * (1.0 - frac) + gamVden * frac
            para[iagamV, ip] = gamVde

            alt_m = alt
            atmos_state = atmos(alt_m, ΔTatmos)
            T0 = atmos_state.T
            p0 = atmos_state.p
            ρ0 = atmos_state.ρ
            a0 = atmos_state.a
            μ0 = atmos_state.μ
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = ρ0
            pare[iemu0, ip] = μ0

            para[iaRange, ip] = R
            para[iaalt, ip] = alt

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
            para[iaWbuoy, ip] = (rhocab - ρ0) * gee * parg[igcabVol]
      end
      para[iaWbuoy, ipdescentn] = 0.0

      # integrate time and weight over descent
      for ip = ipdescent1:ipdescentn

            # velocity calculation from CL, Weight, altitude
            gamVde = para[iagamV, ip]
            cosg = cos(gamVde)
            W = para[iafracW, ip] * WMTO
            BW = W + para[iaWbuoy, ip]
            CL = para[iaCL, ip]
            rho = pare[ierho0, ip]
            V = sqrt(2.0 * BW * cosg / (rho * S * CL))
            Mach = V / pare[iea0, ip]

            para[iaMach, ip] = Mach
            para[iaReunit, ip] = V * rho / pare[iemu0, ip]

            pare[ieu0, ip] = V
            pare[ieM0, ip] = Mach

            # set pitch trim by adjusting CLh
            Wf = W - Wzero
            rfuel = Wf / parg[igWfuel]
            opt_trim_var = TrimVar.CLHtail
            balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

            if (ip == ipdescentn)
                  # use explicitly specified wing cdf,cdp
                  computes_wing_direct = false
            else
                  # use airfoil database for wing cdf,cdp
                  computes_wing_direct = true
            end
            aircraft_drag!(ac, imission, ip, computes_wing_direct)

            # set up for engine calculation
            sing = sin(gamVde)
            cosg = cos(gamVde)
            DoL = para[iaCD, ip] / para[iaCL, ip]
            Fspec = BW * (sing + cosg * DoL)
            pare[ieFe, ip] = Fspec / parg[igneng]

            if initializes_engine
                  pare[iembf, ip] = pare[iembf, ip-1]
                  pare[iemblc, ip] = pare[iemblc, ip-1]
                  pare[iembhc, ip] = pare[iembhc, ip-1]
                  pare[iepif, ip] = pare[iepif, ip-1]
                  pare[iepilc, ip] = pare[iepilc, ip-1]
                  pare[iepihc, ip] = pare[iepihc, ip-1]
                  initializes_engine = false #Apparently, this helps convergence
                                             #on descent, where the engine is at lower throttle

                  # make better estimate for new Tt4, adjusted for new ambient T0
                  dTburn = pare[ieTt4, ip-1] - pare[ieTt3, ip-1]
                  OTR = pare[ieTt3, ip-1] / pare[ieTt2, ip-1]
                  Tt3 = pare[ieT0, ip] * OTR
                  pare[ieTt4, ip] = Tt3 + dTburn + 50.0

                  # make better estimate for new pt5, adjusted for new ambient p0
                  pare[iept5, ip] = pare[iept5, ip-1] * pare[iep0, ip] / pare[iep0, ip-1]

            end

            eng.enginecalc!(ac, "off_design", imission, ip, initializes_engine)

            # store effective thrust, effective TSFC
            F = pare[ieFe, ip] * parg[igneng]
            TSFC = pare[ieTSFC, ip]

            # store integrands for Range and Weight integration
            FoW[ip] = F / (BW * cosg) - DoL

            Vgi[ip] = 1.0 / (V * cosg)

            # if F < 0, then TSFC is not valid, so calculate mdot_fuel directly
            mfuel = pare[iemfuel, ip]
            FFC[ip] = mfuel * gee / (W * V * cosg)

            if (ip > ipdescent1)
                  #  corrector integration step, approximate trapezoidal
                  dh = para[iaalt, ip] - para[iaalt, ip-1]
                  dVsq = pare[ieu0, ip]^2 - pare[ieu0, ip-1]^2

                  FoWavg = 0.5 * (FoW[ip] + FoW[ip-1])
                  FFCavg = 0.5 * (FFC[ip] + FFC[ip-1])
                  Vgiavg = 0.5 * (Vgi[ip] + Vgi[ip-1])

                  dR = para[iaRange, ip] - para[iaRange, ip-1]
                  dt = dR * Vgiavg
                  rW = exp(-dR * FFCavg)

                  para[iatime, ip] = para[iatime, ip-1] + dt
                  para[iafracW, ip] = para[iafracW, ip-1] * rW
            end

            if (ip < ipdescentn)

                  # predictor integration step, forward Euler
                  gamVde = para[iagamV, ip+1]
                  cosg = cos(gamVde)
                  W = para[iafracW, ip+1] * WMTO

                  BW = W + para[iaWbuoy, ip+1]

                  CL = para[iaCL, ip+1]
                  rho = pare[ierho0, ip+1]
                  V = sqrt(2 * BW * cosg / (rho * S * CL))
                  pare[ieu0, ip+1] = V

                  dh = para[iaalt, ip+1] - para[iaalt, ip]
                  dVsq = pare[ieu0, ip+1]^2 - pare[ieu0, ip]^2

                  dR = para[iaRange, ip] - para[iaRange, ip-1]
                  dt = dR * Vgi[ip]
                  rW = exp(-dR * FFC[ip])

                  para[iatime, ip+1] = para[iatime, ip] + dt
                  para[iafracW, ip+1] = para[iafracW, ip] * rW

            end
      end

      # mission fuel fractions and weights
      Wfvent = parm[imWfvent] #Weight of fuel that is vented from tank

      ffvent = Wfvent/WMTO #weight fraction of vented fuel
      
      fracWa = para[iafracW, ipclimb1]
      fracWe = para[iafracW, ipdescentn]
      freserve = parg[igfreserve]
      fburn = fracWa - fracWe + ffvent #include vented fuel, if any
      ffuel = fburn * (1.0 + freserve)
      Wfuel = WMTO * ffuel
      WTO = Wzero + Wfuel

      parm[imWTO] = WTO
      parm[imWfuel] = Wfuel
      #TODO the above calculation does not account for the effect of venting on the flight profile

      # mission PFEI
      Wburn = WMTO * fburn
      parm[imPFEI] = Wburn/gee * parg[igLHVfuel] / (parm[imWpay] * parm[imRange])

end
