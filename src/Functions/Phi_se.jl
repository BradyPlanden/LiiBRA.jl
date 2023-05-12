function Phi_se(Cell, s, z, ElecDef, Φₜ, D, res₀)
    """ 
    Solid-Electrolyte Potential Transfer Function

    Phi_se(Cell, s, z, ElecDef, Φₜ, D, res₀)

    """

    if ElecDef == "Pos"
        Elec = Cell.Pos
    else
        Elec = Cell.Neg
    end

    #Effective Conductivities
    κᵉᶠᶠ = Cell.Const.κ * Elec.ϵ_e^(Elec.κ_brug)
    σᵉᶠᶠ = Elec.σ * Elec.ϵ_s^(Elec.σ_brug)

    #Defining SOC
    θ = Cell.Const.SOC * (Elec.θ_100 - Elec.θ_0) + Elec.θ_0

    #Prepare for j₀
    cs₀ = Elec.cs_max * θ

    #Current Flux Density
    if Cell.Const.CellTyp == "Doyle_94"
        κ = Elec.k_norm / Elec.cs_max / Cell.Const.ce0^(1 - Elec.α)
        j₀ = κ * (Cell.Const.ce0 * (Elec.cs_max - cs₀))^(1 - Elec.α) * cs₀^Elec.α
    else
        j₀ = Elec.k_norm * (Cell.Const.ce0 * cs₀ * (Elec.cs_max - cs₀))^(1 - Elec.α)
    end

    #Resistance + OCP
    Rₜ = R * Cell.Const.T / (j₀ * F^2) + Elec.RFilm
    ∂Uocp = Cell.Const.∂Uocp(ElecDef, θ) / Elec.cs_max

    #Condensing Variable
    ν = @. Elec.L * sqrt((Elec.as / σᵉᶠᶠ + Elec.as / κᵉᶠᶠ) / (Rₜ +
                 ∂Uocp * (Elec.Rs / (F * Elec.Ds)) *
                 (tanh(Elec.β) / (tanh(Elec.β) - Elec.β))))
    ν_∞ = @. Elec.L * sqrt(Elec.as * ((1 / κᵉᶠᶠ) + (1 / σᵉᶠᶠ)) / (Rₜ))

    res₀ .= @. -3 * ∂Uocp / (Elec.as * F * Elec.L * Cell.Const.CC_A * Elec.Rs)

    Φₜ .= @. Elec.L / (Cell.Const.CC_A * ν * sinh(ν)) *
             ((1 / κᵉᶠᶠ) * cosh(ν * z) + (1 / σᵉᶠᶠ) * cosh(ν * (z - 1))) - res₀ / s

    tf₀ = @. (6 * (5 * Elec.Ds * F * Rₜ - ∂Uocp * Elec.Rs) * σᵉᶠᶠ) /
             (30 * Cell.Const.CC_A * Elec.as * Elec.Ds * F * σᵉᶠᶠ * Elec.L) +
             (5 * Elec.as * Elec.Ds * F * Elec.L^2 *
              (σᵉᶠᶠ * (-1 + 3 * z^2) + κᵉᶠᶠ * (2 - 6 * z + 3 * z^2))) /
             (30 * Cell.Const.CC_A * Elec.as * Elec.Ds * F * σᵉᶠᶠ * κᵉᶠᶠ * Elec.L)
    D .= @. Elec.L / (Cell.Const.CC_A * ν_∞ * sinh(ν_∞)) *
            ((1 / κᵉᶠᶠ) * cosh(ν_∞ * z) + (1 / σᵉᶠᶠ) * cosh(ν_∞ * (z - 1))) #D as G(s)->∞
    Φₜ[:, findall(s .== 0)] .= tf₀[:, findall(s .== 0)]
    res₀ .= zeros(length(z))

    if ElecDef == "Pos"
        Φₜ .= -Φₜ
        D .= -D
    end
end
