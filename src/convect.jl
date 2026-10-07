
function convect!(s_interp, s_center, v2_center,  m_center, run_params::RunParams{ConvectionOff, T}) where T
end

# Entropy dips no deeper than this (s is a logarithm) are round-off, not convective zones: conduction
# leaves ~1e-10 wiggles on plateaus flattened at earlier steps, and flattening leaves edge residue at
# the root-finder's precision. They are lifted to the preceding value instead of mixed.
const DS_ROUNDOFF = 1e-9

"Lift round-off-level dips (s[i] below s[i-1] by at most DS_ROUNDOFF) to the preceding value."
function lift_roundoff_dips!(s_center)
    for i in 2:length(s_center)
        d = s_center[i-1] - s_center[i]
        0 < d <= DS_ROUNDOFF && (s_center[i] = s_center[i-1])
    end
end

function convect!(s_interp, s_center, v2_center,  m_center, run_params::RunParams{ConvectionOn, T}) where T
    #while there is a negative entropy gradient

    lift_roundoff_dips!(s_center)
    while !issorted(s_center)
        #find the first region
        println("convection")
        lower_idx, peak_idx, valley_idx, upper_idx = find_convective_zone(s_center)
        # adjust_entropy! root-finds on the interpolant, so it must see the CURRENT s: the passed-in
        # s_interp predates this call, and every earlier pass of this loop has flattened s_center
        s_now = build_s_interp(s_center, log.(m_center))
        adjust_entropy!(s_now, s_center, v2_center, m_center, lower_idx, peak_idx, valley_idx, upper_idx)
        # the flattened zone's edges match the plateau only to root-finder precision; clear any residue
        lift_roundoff_dips!(s_center)
    end
end

function find_convective_zone(s_center)
    n = length(s_center)
    #the geometry is clear in Fig 4 of Stahler, Shu, Taam 1980 ("The hydrostatic core")
    #the idea is that a bracketing range on the true convective zone is from where the uncorrected entropy drops below the local minimum
    #up to where it exceeds the local maximum
    peak_idx = 1 #the index of the first local max of the entropy
    while peak_idx < n && s_center[peak_idx+1] >= s_center[peak_idx]
        peak_idx += 1
    end
    s_peak = s_center[peak_idx] 
    valley_idx = peak_idx
    #now we'll continue iterating to find the local minimum
    while valley_idx < n && s_center[valley_idx+1] <= s_center[valley_idx]
        valley_idx += 1
    end
    s_valley = s_center[valley_idx]

    #Now we find an upper limit on the convective zone, stopping if
    # a. We run out of data
    # b. The entropy exceeds the initial peak 
    upper_idx = valley_idx
    while upper_idx < (n - 1)
        upper_idx += 1

        next_s = s_center[upper_idx]

        # Stop if we exceed the peak entropy
        if next_s > s_peak
            break
        end
        
    end

    #finally, we'll walk towards lower mass from the peak, until s drops below the minimum
    lower_idx = peak_idx
    while lower_idx > 1
        lower_idx -= 1
        next_s = s_center[lower_idx]
        if next_s < s_valley
            break
        end
    end
    return lower_idx, peak_idx, valley_idx, upper_idx
end

function adjust_entropy!(s_interp, s_center, v2_center, m_center, lower_idx, peak_idx, valley_idx, upper_idx)

    logm_valley = log(m_center[valley_idx]) #the mass at the minimum
    logm_peak = log(m_center[peak_idx])
    logm_upper = log(m_center[upper_idx])
    logm_lower = log(m_center[lower_idx])
    s_peak = s_center[peak_idx]
    s_valley = s_center[valley_idx]

    # continuity points m_a (rising branch below the peak) and m_b (rising branch above the valley)
    # where the unmixed profile equals the plateau value s; fall back to the bracket ends
    function zone_edges(s)
        objective(logm) = s_interp(logm) - s
        logm_a = (s_interp(logm_lower) < s < s_interp(logm_peak)) && (lower_idx < peak_idx) ?
            find_zero(objective, (logm_lower, logm_peak)) : logm_lower
        logm_b = (s_interp(logm_upper) > s > s_interp(logm_valley)) && (upper_idx > valley_idx) ?
            find_zero(objective, (logm_valley, logm_upper)) : logm_upper
        return logm_a, logm_b
    end
    # Cell masses (faces are midpoints of the centres; the two end cells are taken symmetric). The
    # energy is measured on exactly the cells that get flattened, so the discrete sum v² s dm is
    # conserved.
    n = length(m_center)
    w = [i == 1 ? m_center[2] - m_center[1] : i == n ? m_center[n] - m_center[n-1] :
         0.5 * (m_center[i+1] - m_center[i-1]) for i in 1:n]
    # Search in log m: a bracket-end fallback is then compared exactly, whereas exp(log(m)) can round
    # past m and drop the end cell (for a zone reaching the core, the peak itself, so convect! never exits)
    logm_center = log.(m_center)
    flattened_cells(s) = begin
        logm_a, logm_b = zone_edges(s)
        searchsortedfirst(logm_center, logm_a):searchsortedlast(logm_center, logm_b)
    end
    # heat defect of flattening those cells to s: sum w v² (s - s_i)  (zero = energy conserved)
    energy_residual(s) = sum(w[i] * v2_center[i] * (s - s_center[i]) for i in flattened_cells(s))

    # The residual is <= 0 at s_valley and >= 0 at s_peak, so bracket the root. (Newton is unsafe:
    # the residual is discontinuous where the zone edges jump to the bracket ends.)
    r_valley = energy_residual(s_valley)
    r_peak = energy_residual(s_peak)
    if !(r_valley <= 0 <= r_peak)
        error("Convective adjustment: energy residual does not bracket a root on [s_valley, s_peak] " *
              "(r_valley = $r_valley, r_peak = $r_peak, zone $lower_idx:$upper_idx)")
    end
    # The residual jumps where a node enters or leaves the flattened set, so the bracketed root can
    # sit on such a jump; once the set is fixed, its plateau is the v²-weighted mean over it, which
    # conserves the discrete invariant exactly (any sub-round-off edge mismatch is lifted by convect!)
    J = flattened_cells(find_zero(energy_residual, (s_valley, s_peak)))
    s_current = sum(w[i] * v2_center[i] * s_center[i] for i in J) / sum(w[i] * v2_center[i] for i in J)

    # Flatten the profile
    s_center[J] .= s_current
    return
end
