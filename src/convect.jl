function convect!(s_interp, s_center, v2_center,  m_center)
    #while there is a negative entropy gradient
    diff(s_center)
    while !issorted(s_center)
        #find the first region 
        
        lower_idx, peak_idx, valley_idx, upper_idx = find_convective_zone(s_center)
        adjust_entropy!(s_interp, s_center, v2_center, m_center, lower_idx, peak_idx, valley_idx, upper_idx)
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
    m_range = (m_center[lower_idx:upper_idx])
    logm_range = log.(m_range)
    v2_range = v2_center[lower_idx:upper_idx]
    s_range = s_center[lower_idx:upper_idx]

    #pre compute the defect integrals
    CDF_v2 = cumul_integrate(m_range, v2_range) #\int_{m_l}^{m_h} v2 dM
    CDF_v2_interp = LinearInterpolation(CDF_v2, logm_range)

    CDF_sv2 = cumul_integrate(m_range, s_range  .* v2_range)#\int_{m_l}^{m_h} s v2 dM
    CDF_sv2_interp = LinearInterpolation(CDF_sv2, logm_range)

    #
    tol = 1e-5
    max_iters=100
    s_current =  s_center[peak_idx] - 1e-4 * s_center[peak_idx]
    logm_a = log(m_center[lower_idx]) #the lower bound on m_a 
    logm_b =log(m_center[upper_idx]) #the upper bound on m_b
    logm_lower = logm_a
    logm_upper = logm_b
    for _ in 1:max_iters
        objective(logm) = s_interp(logm) - s_current
        
        #we first find m_a and m_b by enforcing continuity on s
       
        
        if (s_interp(logm_lower)  < s_current < s_interp(logm_peak))  && (lower_idx < peak_idx) #if we actually have a bracketing interval
            logm_a = find_zero(objective, (logm_lower, logm_peak))
        else

            logm_a = logm_lower
        end

        if (s_interp(logm_upper) > s_current > s_interp(logm_valley)) && (upper_idx > valley_idx)
            logm_b = find_zero(objective, (logm_valley, logm_upper))
        else
            logm_b = logm_upper
        end

        #now we calculate the energy residual
        int_v2_dm = (CDF_v2_interp(logm_b) - CDF_v2_interp(logm_a))
        int_s_v2_dm = (CDF_sv2_interp(logm_b) - CDF_sv2_interp(logm_a))
        residual = s_current * int_v2_dm  - int_s_v2_dm
        #and take a newton step
        delta_s =  residual / int_v2_dm
        s_current -= delta_s
        #and check for convergence
        if residual/int_s_v2_dm < tol
            
            idx_a = searchsortedfirst(m_center, exp(logm_a))
            idx_b = searchsortedlast(m_center, exp(logm_b))
            
            # Flatten the profile
            s_center[idx_a:idx_b] .= s_current
            return 

        end
    end

error("Convective adjustment failed to converge")

end