function Phi_e(Cell, s, z, Φₜ, D, res₀)
    """ 
    Electrolyte Potential Transfer Function

    Phi_e(Cell,s,z,Φₜ,D,res₀)

    """

    #Effective Conductivities
    κᵉᶠᶠₛ = Cell.Const.κ * Cell.Sep.ϵ_e^Cell.Sep.κ_brug
    κᵉᶠᶠₙ = Cell.Const.κ * Cell.Neg.ϵ_e^Cell.Neg.κ_brug
    κᵉᶠᶠₚ = Cell.Const.κ * Cell.Pos.ϵ_e^Cell.Pos.κ_brug
    σᵉᶠᶠₙ = Cell.Neg.σ * Cell.Neg.ϵ_s^Cell.Neg.σ_brug
    σᵉᶠᶠₚ = Cell.Pos.σ * Cell.Pos.ϵ_s^Cell.Pos.σ_brug

    #Defining SOC
    θₙ = Cell.Const.SOC * (Cell.Neg.θ_100 - Cell.Neg.θ_0) + Cell.Neg.θ_0
    θₚ = Cell.Const.SOC * (Cell.Pos.θ_100 - Cell.Pos.θ_0) + Cell.Pos.θ_0

    #Prepare for j₀
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
    νₙ_∞ = @. Cell.Neg.L * sqrt(Cell.Neg.as * ((1 / κᵉᶠᶠₙ) + (1 / σᵉᶠᶠₙ)) / (Rₜⁿ))
    νₚ = @. Cell.Pos.L * sqrt((Cell.Pos.as / σᵉᶠᶠₚ + Cell.Pos.as / κᵉᶠᶠₚ) / (Rₜᵖ +
                  ∂Uocpᵖ * (Cell.Pos.Rs / (F * Cell.Pos.Ds)) *
                  (tanh(Cell.Pos.β) / (tanh(Cell.Pos.β) - Cell.Pos.β))))
    νₚ_∞ = @. Cell.Pos.L * sqrt(Cell.Pos.as * ((1 / κᵉᶠᶠₚ) + (1 / σᵉᶠᶠₚ)) / (Rₜᵖ))

    i = Int64(1)
    # Loop Tf's
    for pt in z
        if pt <= Cell.Neg.L + eps()
            Φₜ[i, :] = @. (Cell.Neg.L * (σᵉᶠᶠₙ / κᵉᶠᶠₙ) *
                           (1 - cosh(νₙ * pt / Cell.Neg.L)) - pt * νₙ * sinh(νₙ) +
                           Cell.Neg.L *
                           (cosh(νₙ) - cosh(νₙ * (Cell.Neg.L - pt) / Cell.Neg.L))) /
                          (Cell.Const.CC_A * (κᵉᶠᶠₙ + σᵉᶠᶠₙ) * sinh(νₙ) * νₙ)
            tf₀ = @. -(pt^2) / (2 * Cell.Const.CC_A * κᵉᶠᶠₙ * Cell.Neg.L)
            D[i] = @. (Cell.Neg.L * (σᵉᶠᶠₙ / κᵉᶠᶠₙ) * (1 - cosh(νₙ_∞ * pt / Cell.Neg.L)) -
                       pt * νₙ_∞ * sinh(νₙ_∞) +
                       Cell.Neg.L *
                       (cosh(νₙ_∞) - cosh(νₙ_∞ * (Cell.Neg.L - pt) / Cell.Neg.L))) /
                      (Cell.Const.CC_A * (κᵉᶠᶠₙ + σᵉᶠᶠₙ) * sinh(νₙ_∞) * νₙ_∞)
            Φₜ[i, findall(s .== 0)] .= tf₀

        elseif pt <= Cell.Neg.L + Cell.Sep.L + eps()
            Φₜ[i, :] = @. (Cell.Neg.L - pt) / (Cell.Const.CC_A * κᵉᶠᶠₛ) +
                          (Cell.Neg.L * ((1 - σᵉᶠᶠₙ / κᵉᶠᶠₙ) * tanh(νₙ / 2) - νₙ)) /
                          (Cell.Const.CC_A * (κᵉᶠᶠₙ + σᵉᶠᶠₙ) * νₙ)
            tf₀ = @. (2 * κᵉᶠᶠₙ * Cell.Neg.L - κᵉᶠᶠₛ * Cell.Neg.L - 2 * κᵉᶠᶠₙ * pt) /
                     (2 * Cell.Const.CC_A * κᵉᶠᶠₙ * κᵉᶠᶠₛ)
            D[i] = @. (Cell.Neg.L - pt) / (Cell.Const.CC_A * κᵉᶠᶠₛ) +
                      (Cell.Neg.L * ((1 - σᵉᶠᶠₙ / κᵉᶠᶠₙ) * tanh(νₙ_∞ / 2) - νₙ_∞)) /
                      (Cell.Const.CC_A * (κᵉᶠᶠₙ + σᵉᶠᶠₙ) * νₙ_∞)
            Φₜ[i, findall(s .== 0)] .= tf₀

        else
            Φₜ[i, :] = @. -Cell.Sep.L / (Cell.Const.CC_A * κᵉᶠᶠₛ) +
                          Cell.Neg.L * ((1 - σᵉᶠᶠₙ / κᵉᶠᶠₙ) * tanh(νₙ / 2) - νₙ) /
                          (Cell.Const.CC_A * (σᵉᶠᶠₙ + κᵉᶠᶠₙ) * νₙ) +
                          (Cell.Pos.L * (-σᵉᶠᶠₚ * cosh(νₚ) +
                            σᵉᶠᶠₚ * cosh((Cell.Const.Ltot - pt) * νₚ / Cell.Pos.L) +
                            κᵉᶠᶠₚ *
                            (cosh((pt - Cell.Neg.L - Cell.Sep.L) * νₚ / Cell.Pos.L) - 1)) -
                           (pt - Cell.Neg.L - Cell.Sep.L) * κᵉᶠᶠₚ * sinh(νₚ) * νₚ) /
                          (Cell.Const.CC_A * κᵉᶠᶠₚ * (κᵉᶠᶠₚ + σᵉᶠᶠₚ) * νₚ * sinh(νₚ))
            tf₀ = @. -Cell.Neg.L / (2 * Cell.Const.CC_A * κᵉᶠᶠₙ) -
                     Cell.Sep.L / (Cell.Const.CC_A * κᵉᶠᶠₛ) -
                     (Cell.Pos.L - (Cell.Const.Ltot - pt)^2 / Cell.Pos.L) /
                     (2 * Cell.Const.CC_A * κᵉᶠᶠₚ)
            D[i] = @. -Cell.Sep.L / (Cell.Const.CC_A * κᵉᶠᶠₛ) +
                      Cell.Neg.L * ((1 - σᵉᶠᶠₙ / κᵉᶠᶠₙ) * tanh(νₙ_∞ / 2) - νₙ_∞) /
                      (Cell.Const.CC_A * (σᵉᶠᶠₙ + κᵉᶠᶠₙ) * νₙ_∞) +
                      (Cell.Pos.L * (-σᵉᶠᶠₚ * cosh(νₚ_∞) +
                        σᵉᶠᶠₚ * cosh((Cell.Const.Ltot - pt) * νₚ_∞ / Cell.Pos.L) +
                        κᵉᶠᶠₚ *
                        (cosh((pt - Cell.Neg.L - Cell.Sep.L) * νₚ_∞ / Cell.Pos.L) - 1)) -
                       (pt - Cell.Neg.L - Cell.Sep.L) * κᵉᶠᶠₚ * sinh(νₚ_∞) * νₚ_∞) /
                      (Cell.Const.CC_A * κᵉᶠᶠₚ * (κᵉᶠᶠₚ + σᵉᶠᶠₚ) * νₚ_∞ * sinh(νₚ_∞))
            Φₜ[i, findall(s .== 0)] .= tf₀
        end
        i += 1
    end

    res₀ .= zeros(length(z))
end
