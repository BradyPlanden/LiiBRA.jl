function Simulate(Cell, Input, Def, Tk, SList, SOC, A₀, B₀, C₀, D₀, t)
    """ 
    Function to simulate generated reduced-order models.

    Sim_Model(Cell, Input, Def, Tk, SList, SOC, A₀, B₀, C₀, D₀, t)

    """
    # Determine time span and allocate arrays
    tlength = size(Input, 1)

    A = Array{Float64}(undef, size(A₀[1]))
    B = Array{Float64}(undef, size(B₀[1]))
    C = Array{Float64}(undef, size(C₀[1]))
    D = Array{Float64}(undef, size(D₀[1]))

    # Selecting SS Models
    υ = findnearest(SList, SOC)
    A = A₀[υ]
    B = B₀[υ]
    C = C₀[υ]
    D = D₀[υ]

    # Capturing Indices
    tfstr = Array{String}(undef, 0, 1)
    for i in 1:length(Cell.Transfer.tfs)
        t1 = t2 = Array{String}(undef, 0, 1)
        for j in 1:length(Cell.Transfer.Locs[i])
            t1 = "$(Cell.Transfer.tfs[i])_$(Cell.Transfer.Elec[i])"
            t2 = [t2; t1]
        end
        tfstr = [tfstr; t2]
    end

    CeInd = findall(isequal("C_e_Na"), tfstr)
    ϕ_ẽInd = findall(isequal("Phi_e_Na"), tfstr)
    CsePosInd = findall(isequal("C_se_Pos"), tfstr)
    CseNegInd = findall(isequal("C_se_Neg"), tfstr)
    ϕ_sNegInd = findall(isequal("Phi_s_Neg"), tfstr)
    ϕ_sPosInd = findall(isequal("Phi_s_Pos"), tfstr)
    ϕ_sePosInd = findall(isequal("Phi_se_Pos"), tfstr)
    ϕ_seNegInd = findall(isequal("Phi_se_Neg"), tfstr)
    FluxNegInd = findall(isequal("Flux_Neg"), tfstr)
    FluxPosInd = findall(isequal("Flux_Pos"), tfstr)

    CeNeg = Cell.Transfer.Locs[1] .<= Cell.Neg.L
    CeNegOffset = Cell.Transfer.Locs[1] .< Cell.Neg.L
    CeSep = Cell.Transfer.Locs[1] .<= Cell.Neg.L + Cell.Sep.L
    CeSepOffset = Cell.Transfer.Locs[1] .< Cell.Neg.L + Cell.Sep.L
    CePos = Cell.Transfer.Locs[1] .<= Cell.Const.Ltot
    CeNegInd = findall(CeNeg .== 1)
    CeSepInd = findall(CeSep .- CeNegOffset .== 1)
    CePosInd = findall(CePos .- CeSepOffset .== 1)

    csegain_neg = C[CseNegInd[1][1], 1] #First Column in C Array (zeros column)
    csegain_pos = C[CsePosInd[1][1], 1] #First Column in C Array (zeros column)

    # Memory Allocation
    Results = (θₙ = Array{Float64}(undef, tlength, 1) .= 0.0,
               θₚ = Array{Float64}(undef, tlength, 1) .= 0.0,
               jeqₙ = Array{Float64}(undef, tlength, 1) .= 0.0,
               jeqₚ = Array{Float64}(undef, tlength, 1) .= 0.0,
               x = Array{Float64}(undef, tlength + 1, size(A, 1)) .= 0.0,
               y = Array{Float64}(undef, tlength, size(C, 1)) .= 0.0,
               Cseₙ = Array{Float64}(undef, tlength, size(CseNegInd, 1)) .= 0.0,
               Cseₚ = Array{Float64}(undef, tlength, size(CsePosInd, 1)) .= 0.0,
               Ce = Array{Float64}(undef, tlength, size(CeInd, 1)) .= Cell.Const.ce0,
               η₀ = Array{Float64}(undef, tlength, 1) .= 0.0,
               ηₙ = Array{Float64}(undef, tlength, size(FluxNegInd, 1)) .= 0.0,
               ηL = Array{Float64}(undef, tlength, 1) .= 0.0,
               ηₚ = Array{Float64}(undef, tlength, size(FluxPosInd, 1)) .= 0.0,
               ϕ_ẽ1 = Array{Float64}(undef, tlength, size(ϕ_ẽInd, 1)) .= 0.0,
               ϕ_ẽ2 = Array{Float64}(undef, tlength, size(CeInd, 1)) .= 0.0,
               ϕse₀ⁿ = Array{Float64}(undef, tlength, 1) .= 0.0, #Replace with length of ϕ_seNegInd @ zero
               jₙ = Array{Float64}(undef, tlength, size(FluxNegInd, 1)) .= 0.0,
               j₀ = Array{Float64}(undef, tlength, 1) .= 0.0,
               jₚ = Array{Float64}(undef, tlength, size(FluxPosInd, 1)) .= 0.0,
               jL = Array{Float64}(undef, tlength, 1) .= 0.0,
               Rₜⁿ = Array{Float64}(undef, tlength, 1) .= 0.0,
               Rₜᵖ = Array{Float64}(undef, tlength, 1) .= 0.0,
               Uocpⁿ = Array{Float64}(undef, tlength, 1) .= 0.0,
               Uocpᵖ = Array{Float64}(undef, tlength, 1) .= 0.0,
               Cell_V = Array{Float64}(undef, tlength, 1) .= 0.0,
               ϕ_e = Array{Float64}(undef, tlength, size(CeInd, 1)) .= 0.0,
               Cell_SOC = Array{Float64}(undef, tlength, 1) .= 0,
               Iapp = Array{Float64}(undef, tlength + 1, 1) .= 0,
               t = t,
               tₑ = Int64(1))

    # Defining SOC
    SOCₙ = SOC * (Cell.Neg.θ_100 - Cell.Neg.θ_0) + Cell.Neg.θ_0
    SOCₚ = SOC * (Cell.Pos.θ_100 - Cell.Pos.θ_0) + Cell.Pos.θ_0
    Results.θₙ[1] = SOCₙ
    Results.θₚ[1] = SOCₚ
    tₑ = Int64(1)

    # Loop through time - compute dependent variables (voltage, flux, etc.) #
    for i in 0:(tlength - 1)
        cs_neg_avg = Results.x[i + 1, 1] * csegain_neg + SOCₙ * Cell.Neg.cs_max < 0.0 ?
                     0.0 : Results.x[i + 1, 1] * csegain_neg + SOCₙ * Cell.Neg.cs_max #Zero if < 0
        cs_pos_avg = Results.x[i + 1, 1] * csegain_pos + SOCₚ * Cell.Pos.cs_max < 0.0 ?
                     0.0 : Results.x[i + 1, 1] * csegain_pos + SOCₚ * Cell.Pos.cs_max #Zero if < 0

        if cs_neg_avg > Cell.Neg.cs_max
            cs_neg_avg = Cell.Neg.cs_max
        end

        if cs_pos_avg > Cell.Pos.cs_max
            cs_pos_avg = Cell.Pos.cs_max
        end

        # Reaction Rates
        if Cell.Const.CellTyp == "Doyle_94"
            kₙ = Cell.Neg.k_norm
            kₚ = Cell.Pos.k_norm
            Results.jeqₙ[i + 1] = kₙ * sqrt(cs_neg_avg * Cell.Const.ce0 *
                                       (Cell.Neg.cs_max - cs_neg_avg))
            Results.jeqₚ[i + 1] = kₚ * sqrt(cs_pos_avg * Cell.Const.ce0 *
                                       (Cell.Pos.cs_max - cs_pos_avg))
        else
            kₙ = Cell.Neg.k_norm
            kₚ = Cell.Pos.k_norm
            Results.jeqₙ[i + 1] = Cell.Neg.k_norm *
                                  (Cell.Const.ce0 * cs_neg_avg *
                                   (Cell.Neg.cs_max - cs_neg_avg))^(1 - Cell.Neg.α)
            Results.jeqₚ[i + 1] = Cell.Pos.k_norm *
                                  (Cell.Const.ce0 * cs_pos_avg *
                                   (Cell.Pos.cs_max - cs_pos_avg))^(1 - Cell.Pos.α)
        end

        Results.θₙ[i + 1] = cs_neg_avg / Cell.Neg.cs_max
        Results.θₚ[i + 1] = cs_pos_avg / Cell.Pos.cs_max
        Results.Cell_SOC[i + 1] = (Results.θₙ[i + 1] - Cell.Neg.θ_0) /
                                  (Cell.Neg.θ_100 - Cell.Neg.θ_0)

        javg_neg = Results.Iapp[i + 1] / (Cell.Neg.as * F * Cell.Neg.L * Cell.Const.CC_A)
        javg_pos = Results.Iapp[i + 1] / (Cell.Pos.as * F * Cell.Pos.L * Cell.Const.CC_A)

        Arr_Factor = ((1 / Cell.Const.T_ref) - (1 / Tk[i + 1])) / R
        κneg = Cell.Const.κf(mean(Results.Ce[i + 1, CeNegInd])) *
               exp(Cell.Const.Ea_κ * Arr_Factor)
        κpos = Cell.Const.κf(mean(Results.Ce[i + 1, CePosInd])) *
               exp(Cell.Const.Ea_κ * Arr_Factor)
        κsep = Cell.Const.κf(mean(Results.Ce[i + 1, CeSepInd])) *
               exp(Cell.Const.Ea_κ * Arr_Factor)
        σ_neg = Cell.Neg.σ * exp(Cell.Const.Ea_κ * Arr_Factor)
        σ_pos = Cell.Pos.σ * exp(Cell.Const.Ea_κ * Arr_Factor)
        κ_eff_Neg = κneg * (Cell.Neg.ϵ_e^(Cell.Neg.κ_brug))
        κ_eff_Sep = κsep * (Cell.Sep.ϵ_e^(Cell.Sep.κ_brug))
        κ_eff_Pos = κpos * (Cell.Pos.ϵ_e^(Cell.Pos.κ_brug))
        σ_eff_Neg = σ_neg * Cell.Neg.ϵ_s^Cell.Neg.σ_brug
        σ_eff_Pos = σ_pos * Cell.Pos.ϵ_s^Cell.Pos.σ_brug

        # Resistances
        Results.Rₜⁿ[i + 1] = (Tk[i + 1] * R) /
                             (F^2 * sqrt(Results.jeqₙ[i + 1]^2 + javg_neg^2 / 4)) +
                             Cell.Neg.RFilm
        Results.Rₜᵖ[i + 1] = (Tk[i + 1] * R) /
                             (F^2 * sqrt(Results.jeqₚ[i + 1]^2 + javg_pos^2 / 4)) +
                             Cell.Pos.RFilm

        # Condensing Variable
        ν_neg = @. Cell.Neg.L * sqrt((Cell.Neg.as * (1 / κ_eff_Neg + 1 / σ_eff_Neg)) /
                        Results.Rₜⁿ[i + 1])
        ν_pos = @. Cell.Pos.L * sqrt((Cell.Pos.as * (1 / κ_eff_Pos + 1 / σ_eff_Pos)) /
                        Results.Rₜᵖ[i + 1])

        # Relinearise dependent on ν, σ, κ
        D = D_Linear(Cell, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos,
                     κ_eff_Sep)

        # Interpolate C & D Matrices
        C = interp(C₀, SList, Results.Cell_SOC[i + 1])
        # D = interp(D₀,SList,Results.Cell_SOC[i+1])

        # SS Output
        Results.y[i + 1, :] = C * Results.x[i + 1, :] + D * Results.Iapp[i + 1]

        # Concentrations & Force electrode concentration maximum
        Results.Cseₙ[i + 1, :] = (SOCₙ .* Cell.Neg.cs_max .+
                                  Results.y[i + 1, CseNegInd]) >
                                 ones(size(Results.Cseₙ, 2)) * Cell.Neg.cs_max ?
                                 ones(size(Results.Cseₙ, 2)) * Cell.Neg.cs_max :
                                 (SOCₙ .* Cell.Neg.cs_max .+
                                  Results.y[i + 1, CseNegInd])
        Results.Cseₚ[i + 1, :] = (SOCₚ .* Cell.Pos.cs_max .+
                                  Results.y[i + 1, CsePosInd]) >
                                 ones(size(Results.Cseₚ, 2)) * Cell.Pos.cs_max ?
                                 ones(size(Results.Cseₚ, 2)) * Cell.Pos.cs_max :
                                 (SOCₚ .* Cell.Pos.cs_max .+
                                  Results.y[i + 1, CsePosInd])
        Results.Ce[i + 1, :] = @. Cell.Const.ce0 + Results.y[i + 1, CeInd]

        # Potentials
        Results.Uocpⁿ[i + 1] = Cell.Const.Uocp("Neg",
                                               Results.Cseₙ[i + 1, 1] /
                                               Cell.Neg.cs_max)
        Results.Uocpᵖ[i + 1] = Cell.Const.Uocp("Pos",
                                               Results.Cseₚ[i + 1, 1] /
                                               Cell.Pos.cs_max)
        Results.ϕse₀ⁿ[i + 1] = Results.y[i + 1, ϕ_seNegInd[1]] + Results.Uocpⁿ[i + 1] #Location 0
        Results.ϕ_ẽ1[i + 1, :] = Results.y[i + 1, ϕ_ẽInd]

        Results.ϕ_ẽ2[i + 1, :] = @. ((Tk[i + 1] * 2 * R *
                                       (1 - Cell.Const.tpf(Results.Ce[i + 1, :]))) / F) *
                                     (log(Results.Ce[i + 1, :] / Results.Ce[i + 1, 1]))
        Results.ϕ_e[i + 1, :] = @. [0; Results.ϕ_ẽ1[i + 1, :]] + Results.ϕ_ẽ2[i + 1, :] -
                                   Results.ϕse₀ⁿ[i + 1]

        # Flux
        Results.jₙ[i + 1, :] = Results.y[i + 1, FluxNegInd]
        Results.j₀[i + 1] = Results.y[i + 1, FluxNegInd[1]]
        Results.jₚ[i + 1, :] = Results.y[i + 1, FluxPosInd]
        Results.jL[i + 1] = Results.y[i + 1, FluxPosInd[1]]

        # Neg
        j₀ⁿ = findmax([ones(size(Results.Cseₙ, 2)) * eps() (kₙ .*
                                                            (Results.Cseₙ[i + 1,
                                                                          :] .^
                                                             Cell.Neg.α) .*
                                                            (Results.Ce[i + 1, 1] .^
                                                             (1 - Cell.Neg.α))) .*
                                                           (Cell.Neg.cs_max .-
                                                            Results.Cseₙ[i + 1, :]) .^
                                                           (1 - Cell.Neg.α)], dims = 2)[1]
        Results.η₀[i + 1] = Tk[i + 1] * 2 * R / F *
                            asinh(Results.j₀[i + 1] / (2 * j₀ⁿ[1]))
        Results.ηₙ[i + 1, :] = @. (Tk[i + 1] * 2 * R) / F *
                                  asinh((Results.jₙ[i + 1, :]) / (2 * j₀ⁿ))

        # Pos
        j₀ᵖ = findmax([ones(size(Results.Cseₚ, 2)) * eps() (kₚ .*
                                                            (Results.Cseₚ[i + 1,
                                                                          :] .^
                                                             Cell.Pos.α) .*
                                                            (Results.Ce[i + 1, 1] .^
                                                             (1 - Cell.Pos.α))) .*
                                                           (Cell.Pos.cs_max .-
                                                            Results.Cseₚ[i + 1, :]) .^
                                                           (1 - Cell.Pos.α)], dims = 2)[1]
        Results.ηL[i + 1] = (Tk[i + 1] * 2 * R) / F *
                            asinh(Results.jL[i + 1] / (2 * j₀ᵖ[1]))
        Results.ηₚ[i + 1, :] = @. (Tk[i + 1] * 2 * R) / F *
                                  asinh(Results.jₚ[i + 1, :] / (2 * j₀ᵖ))

        # Cell Voltage
        Results.Cell_V[i + 1] = @. (Results.Uocpᵖ[i + 1] - Results.Uocpⁿ[i + 1]) +
                                   (Results.ηL[i + 1] - Results.η₀[i + 1]) +
                                   (Results.ϕ_ẽ1[i + 1, end] + Results.ϕ_ẽ2[i + 1, end]) +
                                   (Cell.Pos.RFilm * Results.jL[i + 1] -
                                    Cell.Neg.RFilm * Results.j₀[i + 1]) * F

        # ϕ_s
        ϕ_s_neg = Results.y[i + 1, ϕ_sNegInd]
        ϕ_s_pos = @. Results.y[i + 1, ϕ_sPosInd] + Results.Cell_V[i + 1]

        # Interpolate A Matrix
        A = interp(A₀, SList, Results.Cell_SOC[i + 1])
        B = interp(B₀, SList, Results.Cell_SOC[i + 1])

        # Update States
        Results.x[i + 2, :] = A * Results.x[i + 1, :] + B * Results.Iapp[i + 1]

        if Def == "Power"
            Results.Iapp[i + 2] = Input[i + 1, 1] / Results.Cell_V[i + 1]
        else
            Results.Iapp[i + 2] = Input[i + 1, 1]
        end

        if Results.Cell_V[i + 1] <= Cell.Const.Vmin
            break
        end

        tₑ = i
    end

    return Results, tₑ
end
