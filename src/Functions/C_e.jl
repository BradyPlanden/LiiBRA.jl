function C_e(Cell, s, z, tf, D, res₀)
    """ 
    Electrolyte Concentration Transfer Function

    C_e(Cell, s, z, tf, D, res₀)

    """

    #Effective Conductivities
    ζ = (1 - Cell.Const.t_plus) / F
    κᵉᶠᶠₙ = Cell.Const.κ * Cell.Neg.ϵ_e^Cell.Neg.κ_brug
    κᵉᶠᶠₚ = Cell.Const.κ * Cell.Pos.ϵ_e^Cell.Pos.κ_brug
    σᵉᶠᶠₙ = Cell.Neg.σ * Cell.Neg.ϵ_s^Cell.Neg.σ_brug
    σᵉᶠᶠₚ = Cell.Pos.σ * Cell.Pos.ϵ_s^Cell.Pos.σ_brug

    #Defining SOC
    θₙ = Cell.Const.SOC * (Cell.Neg.θ_100 - Cell.Neg.θ_0) + Cell.Neg.θ_0
    θₚ = Cell.Const.SOC * (Cell.Pos.θ_100 - Cell.Pos.θ_0) + Cell.Pos.θ_0

    #Prepare for j0
    cs₀ⁿ = Cell.Neg.cs_max * θₙ
    cs₀ᵖ = Cell.Pos.cs_max * θₚ

    #Current Flux Density
    if Cell.Const.CellTyp == "Doyle_94"
        κₚ = Cell.Pos.k_norm / Cell.Pos.cs_max / Cell.Const.ce0^(1 - Cell.Pos.α)
        κₙ = Cell.Neg.k_norm / Cell.Neg.cs_max / Cell.Const.ce0^(1 - Cell.Neg.α)
        j₀ⁿ = κₙ * (Cell.Const.ce0 * (Cell.Neg.cs_max - cs₀ⁿ))^(1 - Cell.Neg.α) *
              cs₀ⁿ^Cell.Neg.α
        j₀ᵖ = κₚ * (Cell.Const.ce0 * (Cell.Pos.cs_max - cs₀ᵖ))^(1 - Cell.Pos.α) *
              cs₀ᵖ^Cell.Pos.α
    else
        j₀ⁿ = Cell.Neg.k_norm *
              (Cell.Const.ce0 * cs₀ⁿ * (Cell.Neg.cs_max - cs₀ⁿ))^(1 - Cell.Neg.α)
        j₀ᵖ = Cell.Pos.k_norm *
              (Cell.Const.ce0 * cs₀ᵖ * (Cell.Pos.cs_max - cs₀ᵖ))^(1 - Cell.Pos.α)
    end

    #Resistances
    Rₜⁿ = R * Cell.Const.T / (j₀ⁿ * F^2) + Cell.Neg.RFilm
    Rₜᵖ = R * Cell.Const.T / (j₀ᵖ * F^2) + Cell.Pos.RFilm

    #OCP derivative
    ∂Uocpⁿ = Cell.Const.∂Uocp("Neg", θₙ) / Cell.Neg.cs_max
    ∂Uocpᵖ = Cell.Const.∂Uocp("Pos", θₚ) / Cell.Pos.cs_max

    #Condensing Variable
    νₙ = @. Cell.Neg.L * sqrt((Cell.Neg.as / σᵉᶠᶠₙ + Cell.Neg.as / κᵉᶠᶠₙ) / (Rₜⁿ +
                  ∂Uocpⁿ * (Cell.Neg.Rs / (F * Cell.Neg.Ds)) *
                  (tanh(Cell.Neg.β) / (tanh(Cell.Neg.β) - Cell.Neg.β))))

    νₚ = @. Cell.Pos.L * sqrt((Cell.Pos.as / σᵉᶠᶠₚ + Cell.Pos.as / κᵉᶠᶠₚ) / (Rₜᵖ +
                  ∂Uocpᵖ * (Cell.Pos.Rs / (F * Cell.Pos.Ds)) *
                  (tanh(Cell.Pos.β) / (tanh(Cell.Pos.β) - Cell.Pos.β))))

    R_ce = find_zeros(x -> flambda(Cell, x), float(0), Cell.Const.CeRootRange)
    if size(R_ce, 1) >= Cell.Const.Ce_M + 1
        λ = R_ce[2:(Cell.Const.Ce_M + 1)]
    else
        throw(DomainError(R_ce,
                          "Ce roots 'R_ce' doesn't contain enough values, increase Cell.Const.CeRootRange"))
    end

    #Create all k's
    in₁ = @. sqrt(λ * Cell.Neg.ϵ_e / Cell.Const.D1)
    in₂ = @. sqrt(λ * Cell.Sep.ϵ_e / Cell.Const.D2)
    in₃ = @. sqrt(λ * Cell.Pos.ϵ_e / Cell.Const.D3)

    Bndⁿ₁ = Cell.Neg.L * in₁
    Bndˢ₀ = Cell.Neg.L * in₂
    Bndˢ₁ = (Cell.Neg.L + Cell.Sep.L) * in₂
    Bndᵖ₀ = (Cell.Neg.L + Cell.Sep.L) * in₃
    Bndᵖ₁ = (Cell.Neg.L + Cell.Sep.L + Cell.Pos.L) * in₃
    Bndᵖ₂ = Cell.Pos.L * in₃

    #Scaled coefficients
    k3ₛ = @. cos(Bndⁿ₁) * cos(Bndˢ₀) +
             Cell.Const.D1 * in₁ * sin(Bndⁿ₁) * sin(Bndˢ₀) / (Cell.Const.D2 * in₂)
    k4ₛ = @. cos(Bndⁿ₁) * sin(Bndˢ₀) -
             Cell.Const.D1 * in₁ * cos(Bndˢ₀) * sin(Bndⁿ₁) / (Cell.Const.D2 * in₂)
    k5ₛ = @. k3ₛ * (cos(Bndˢ₁) * cos(Bndᵖ₀) +
              Cell.Const.D2 * in₂ * sin(Bndˢ₁) * sin(Bndᵖ₀) / (Cell.Const.D3 * in₃)) +
             k4ₛ * (sin(Bndˢ₁) * cos(Bndᵖ₀) -
              Cell.Const.D2 * in₂ * cos(Bndˢ₁) * sin(Bndᵖ₀) / (Cell.Const.D3 * in₃))
    k6ₛ = @. k3ₛ * (cos(Bndˢ₁) * sin(Bndᵖ₀) -
              Cell.Const.D2 * in₂ * sin(Bndˢ₁) * cos(Bndᵖ₀) / (Cell.Const.D3 * in₃)) +
             k4ₛ * (sin(Bndˢ₁) * sin(Bndᵖ₀) +
              Cell.Const.D2 * in₂ * cos(Bndˢ₁) * cos(Bndᵖ₀) /
              (Cell.Const.D3 * in₃))

    #Solving for k1:
    ψⁱⁿᵗ₁ = @. Cell.Neg.ϵ_e * (2 * Bndⁿ₁ + sin(2 * Bndⁿ₁)) / (4 * in₁)
    ψⁱⁿᵗ₂ = @. Cell.Sep.ϵ_e / (4 * in₂) * (2 * (k3ₛ^2 + k4ₛ^2) * Cell.Sep.L * in₂ +
                2 * k3ₛ * k4ₛ * cos(2 * Bndˢ₀) -
                2 * k3ₛ * k4ₛ * cos(2 * Bndˢ₁) -
                (k3ₛ - k4ₛ) * (k3ₛ + k4ₛ) *
                (sin(2 * Bndˢ₀) - sin(2 * Bndˢ₁)))
    ψⁱⁿᵗ₃ = @. Cell.Pos.ϵ_e / (4 * in₃) * (2 * (k5ₛ^2 + k6ₛ^2) * Cell.Pos.L * in₃ +
                2 * k5ₛ * k6ₛ * cos(2 * Bndᵖ₀) -
                2 * k5ₛ * k6ₛ * cos(2 * Bndᵖ₁) -
                (k5ₛ - k6ₛ) * (k5ₛ + k6ₛ) *
                (sin(2 * Bndᵖ₀) - sin(2 * Bndᵖ₁)))

    k1 = @. 1 / sqrt(ψⁱⁿᵗ₁ + ψⁱⁿᵗ₂ + ψⁱⁿᵗ₃)
    k3 = @. k1 * k3ₛ
    k4 = @. k1 * k4ₛ
    k5 = @. k1 * k5ₛ
    k6 = @. k1 * k6ₛ

    jₙ = @. k1 * ζ * νₙ *
            (Bndⁿ₁ * (κᵉᶠᶠₙ + σᵉᶠᶠₙ * cosh(νₙ)) * sin(Bndⁿ₁) +
             (κᵉᶠᶠₙ + σᵉᶠᶠₙ * cos(Bndⁿ₁)) * sinh(νₙ) * νₙ) /
            (Cell.Const.CC_A * (κᵉᶠᶠₙ + σᵉᶠᶠₙ) * (Bndⁿ₁^2 + νₙ^2) * sinh(νₙ))
    tf₀ⁿ = @. k1 * ζ * sin(Bndⁿ₁) / (Cell.Const.CC_A * Bndⁿ₁)
    jₙ[:, findall(s .== 0)] .= tf₀ⁿ[:, findall(s .== 0)]

    jₚ = @. -ζ * νₚ / (Cell.Const.CC_A * (κᵉᶠᶠₚ + σᵉᶠᶠₚ) * (Bndᵖ₂^2 + νₚ^2) *
                       sinh(νₚ)) * (-k6 * Bndᵖ₂ * cos(Bndᵖ₁) * (σᵉᶠᶠₚ + κᵉᶠᶠₚ * cosh(νₚ)) +
             Bndᵖ₂ * (κᵉᶠᶠₚ + σᵉᶠᶠₚ * cosh(νₚ)) * (k6 * cos(Bndᵖ₀) - k5 * sin(Bndᵖ₀)) +
             k5 * Bndᵖ₂ * (σᵉᶠᶠₚ + κᵉᶠᶠₚ * cosh(νₚ)) * sin(Bndᵖ₁) +
             sinh(νₚ) *
             (k5 * σᵉᶠᶠₚ * cos(Bndᵖ₀) + k5 * κᵉᶠᶠₚ * cos(Bndᵖ₁) +
              k6 * σᵉᶠᶠₚ * sin(Bndᵖ₀) + k6 * κᵉᶠᶠₚ * sin(Bndᵖ₁)) * νₚ)
    tf₀ᵖ = @. -ζ * (k6 * (cos(Bndᵖ₀) - cos(Bndᵖ₁)) +
                    k5 * (sin(Bndᵖ₁) - sin(Bndᵖ₀))) /
              (Cell.Const.CC_A * Bndᵖ₂)
    jₚ[:, findall(s .== 0)] .= tf₀ᵖ[:, findall(s .== 0)]

    ψ = fill(float(0), length(z), length(λ))
    for lp in 1:length(λ)
        i = Int64(1)
        for x in z #Eigen Weighting
            if x < Cell.Neg.L + eps()
                ψ[i, lp] = k1[lp] * cos(in₁[lp] * x) #negative
            elseif x > (Cell.Neg.L + Cell.Sep.L) - eps()
                ψ[i, lp] = k5[lp] * cos(in₃[lp] * x) + k6[lp] * sin(in₃[lp] * x) # postive
            else
                ψ[i, lp] = k3[lp] * cos(in₂[lp] * x) + k4[lp] * sin(in₂[lp] * x) # separator
            end
            i = i + 1
        end
    end

    tf .= ψ * ((jₙ .+ jₚ) ./ (s .+ λ))
    D .= zeros(length(z))
    res₀ .= zeros(length(z))
end

function flambda(Cell, λ)
    Lⁿₛ = Cell.Const.Lnegsep
    s₁ = sqrt(λ * Cell.Neg.ϵ_e / Cell.Const.D1)
    s₂ = sqrt(λ * Cell.Sep.ϵ_e / Cell.Const.D2)
    s₃ = sqrt(λ * Cell.Pos.ϵ_e / Cell.Const.D3)
    k3 = @. (cos(s₁ * Cell.Neg.L) * cos(s₂ * Cell.Neg.L) +
             Cell.Const.D1 * s₁ * sin(s₁ * Cell.Neg.L) * sin(s₂ * Cell.Neg.L) /
             (Cell.Const.D2 * s₂))
    k4 = @. (cos(s₁ * Cell.Neg.L) * sin(s₂ * Cell.Neg.L) -
             Cell.Const.D1 * s₁ * cos(s₂ * Cell.Neg.L) * sin(s₁ * Cell.Neg.L) /
             (Cell.Const.D2 * s₂))
    k5 = @. k3 * (cos(s₂ * Lⁿₛ) * cos(s₃ * Lⁿₛ) +
             Cell.Const.D2 * s₂ * sin(s₂ * Lⁿₛ) * sin(s₃ * Lⁿₛ) /
             (Cell.Const.D3 * s₃)) +
            k4 * (sin(s₂ * Lⁿₛ) * cos(s₃ * Lⁿₛ) -
             Cell.Const.D2 * s₂ * cos(s₂ * Lⁿₛ) * sin(s₃ * Lⁿₛ) /
             (Cell.Const.D3 * s₃))
    k6 = @. k3 * (cos(s₂ * Lⁿₛ) * sin(s₃ * Lⁿₛ) -
             Cell.Const.D2 * s₂ * sin(s₂ * Lⁿₛ) * cos(s₃ * Lⁿₛ) /
             (Cell.Const.D3 * s₃)) +
            k4 * (sin(s₂ * Lⁿₛ) * sin(s₃ * Lⁿₛ) +
             Cell.Const.D2 * s₂ * cos(s₂ * Lⁿₛ) * cos(s₃ * Lⁿₛ) /
             (Cell.Const.D3 * s₃))

    if λ == 0
        return zero(typeof(λ))
    else
        return @. -k5 * s₃ * sin(s₃ * Cell.Const.Ltot) +
                  k6 * s₃ * cos(s₃ * Cell.Const.Ltot)
    end
end
