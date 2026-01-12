
@kwdef struct RunParams{QuantitySigma, QuantityL, QuantityV, QuantityT, QuantityRho}
    C::Float64 # Newton's Constant (G) in compatible units (float magnitude)
    σ0::QuantitySigma
    w::QuantityV
    n::Float64
    t_end::QuantityT
    rho_initial::FunctionWrapper{QuantityRho, Tuple{QuantityL}} 
    r_max::QuantityL
    r_min::QuantityL
end