function C_se(Cell, s, z, ElecDef, cseₜ, D, res₀)
    """ 
    Concentration Solid-Electrolyte Transfer Function

    C_se(Cell,s,z,ElecDef)

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

    #Resistance & OCP
    Rₜ = R * Cell.Const.T / (j₀ * F^2) + Elec.RFilm
    ∂Uocp = Cell.Const.∂Uocp(ElecDef, θ) / Elec.cs_max

    #Condensing Variable
    ν = @. Elec.L * sqrt((Elec.as / σᵉᶠᶠ + Elec.as / κᵉᶠᶠ) / (Rₜ +
                 ∂Uocp * (Elec.Rs / (F * Elec.Ds)) *
                 (tanh(Elec.β) / (tanh(Elec.β) - Elec.β))))

    res₀ .= -3 / (Elec.as * F * Elec.L * Cell.Const.CC_A * Elec.Rs) *
            ones(length(z)) #Residual Variable for Pole Removal

    cseₜ .= @. ν * Elec.Rs * (σᵉᶠᶠ * cosh(ν * z) + κᵉᶠᶠ * cosh(ν * (z - 1))) *
               tanh(Elec.β) /
               (Elec.as * F * Elec.L * Cell.Const.CC_A * Elec.Ds *
                (κᵉᶠᶠ + σᵉᶠᶠ) * sinh(ν) * (tanh(Elec.β) - Elec.β)) - res₀ / s

    tf₀ = @. (5 * Elec.as * Elec.Ds * F * Elec.L^2 *
              (κᵉᶠᶠ * (2 - 6 * z + 3 * z^2) + (3 * z^2 - 1) * σᵉᶠᶠ) -
              6 * ∂Uocp * Elec.Rs * κᵉᶠᶠ * σᵉᶠᶠ) /
             (30 * Cell.Const.CC_A * Elec.as * Elec.Ds * ∂Uocp * F *
              Elec.L * κᵉᶠᶠ * σᵉᶠᶠ)

    cseₜ[:, findall(s .== 0)] .= tf₀[:, findall(s .== 0)]
    D .= zeros(length(z))

    if ElecDef == "Pos"
        cseₜ .= -cseₜ
        res₀ .= -res₀
    end

    if abs.(cseₜ[:, 1]) > abs.(cseₜ[:, 2]) * 10
        cseₜ[:, 1] = cseₜ[:, 2] * 10
    end
end
