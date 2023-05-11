function Phi_s(Cell, s, z, ElecDef, Φₜ, D, res₀)
    """ 
    Solid Potential Transfer Function

    Phi_s(Cell, s, z, ElecDef, Φₜ, D, res₀)

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
    ν_∞ = @. Elec.L * sqrt((Elec.as / κᵉᶠᶠ + Elec.as / σᵉᶠᶠ) / (Rₜ))

    Φₜ .= @. (-Elec.L * (κᵉᶠᶠ * (cosh(ν) - cosh((z - 1) * ν))) -
              Elec.L * (σᵉᶠᶠ * (1 - cosh(z * ν) + z * ν * sinh(ν)))) /
             (Cell.Const.CC_A * σᵉᶠᶠ * (κᵉᶠᶠ + σᵉᶠᶠ) * ν * sinh(ν))
    D .= @. (-Elec.L * (κᵉᶠᶠ * (cosh(ν_∞) - cosh((z - 1) * ν_∞))) -
             Elec.L * (σᵉᶠᶠ * (1 - cosh(z * ν_∞) + z * ν_∞ * sinh(ν_∞)))) /
            (Cell.Const.CC_A * σᵉᶠᶠ * (κᵉᶠᶠ + σᵉᶠᶠ) * ν_∞ * sinh(ν_∞)) #D as G->∞
    tf₀ = @. Elec.L * (z - 2) * z / (2 * Cell.Const.CC_A * σᵉᶠᶠ)
    Φₜ[:, findall(s .== 0)] .= tf₀[:, findall(s .== 0)]
    res₀ .= zeros(length(z))

    if ElecDef == "Pos"
        Φₜ .= -Φₜ
        D .= -D
    end
end
