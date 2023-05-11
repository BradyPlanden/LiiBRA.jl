function Flux(Cell, s, z, ElecDef, jₜ, D, res₀)
    """ 
    Flux Transfer Function

    Flux(Cell,s,z,ElecDef, jₜ, D, res₀)

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
    ν̂ = @. Elec.L * sqrt(Elec.as * ((1 / κᵉᶠᶠ) + (1 / σᵉᶠᶠ)) / (Rₜ))

    #Transfer Function
    jₜ .= @. ν * (σᵉᶠᶠ * cosh(ν * z) + κᵉᶠᶠ * cosh(ν * (z - 1))) /
             (Elec.as * F * Elec.L * Cell.Const.CC_A * (κᵉᶠᶠ + σᵉᶠᶠ) * sinh(ν))

    D .= @. ν̂ * (σᵉᶠᶠ * cosh(ν̂ * z) + κᵉᶠᶠ * cosh(ν̂ * (z - 1))) /
            (Elec.as * F * Elec.L * Cell.Const.CC_A * (κᵉᶠᶠ + σᵉᶠᶠ) * sinh(ν̂))

    tf₀ = ones(size(z, 1)) * 1 / (Cell.Const.CC_A * Elec.as * F * Elec.L)
    jₜ[:, findall(s .== 0)] .= tf₀[:, findall(s .== 0)]
    res₀ .= zeros(length(z))

    if ElecDef == "Pos"
        jₜ .= -jₜ
        D .= -D
    end
end
