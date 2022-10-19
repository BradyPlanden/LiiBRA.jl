function Sim_Model(Cell,Input,Def,Tk,SList,SOC,A0,B0,C0,D0,t)
    """ 
    Function to simulate generated reduced-order models.

    Sim_Model(Cell,Input,Def,Tk,SOC,A,B,C,D)
    
    """
    # Determine time span and allocate arrays
    tlength = size(Input,1)

    A = Array{Float64}(undef,size(A0[1]))
    B = Array{Float64}(undef,size(B0[1]))
    C = Array{Float64}(undef,size(C0[1]))
    D = Array{Float64}(undef,size(D0[1]))

    # Selecting SS Models
    υ = findnearest(SList,SOC)
    A = A0[υ]
    B = B0[υ]
    C = C0[υ]
    D = D0[υ]

    # Capturing Indices
    tfstr = Array{String}(undef,0,1)
    for i in 1:length(Cell.Transfer.tfs)
        t1 = t2 = Array{String}(undef,0,1)
        for j in 1:length(Cell.Transfer.Locs[i])
            t1 ="$(Cell.Transfer.tfs[i])_$(Cell.Transfer.Elec[i])"
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


    CeNeg = Cell.Transfer.Locs[1].<=Cell.Neg.L
    CeNegOffset = Cell.Transfer.Locs[1].<Cell.Neg.L
    CeSep = Cell.Transfer.Locs[1].<=Cell.Neg.L+Cell.Sep.L
    CeSepOffset = Cell.Transfer.Locs[1].<Cell.Neg.L+Cell.Sep.L
    CePos = Cell.Transfer.Locs[1].<=Cell.Const.Ltot
    CeNegInd = findall(CeNeg .== 1)
    CeSepInd = findall(CeSep.-CeNegOffset .== 1)
    CePosInd = findall(CePos.-CeSepOffset .== 1)

    csegain_neg = C[CseNegInd[1][1],1] #First Column in C Array (zeros column)
    csegain_pos = C[CsePosInd[1][1],1] #First Column in C Array (zeros column)

    # Memory Allocation
    Results =   (
                θ_neg = Array{Float64}(undef,tlength,1) .= 0.,
                θ_pos = Array{Float64}(undef,tlength,1) .= 0.,
                jeq_neg = Array{Float64}(undef,tlength,1) .= 0.,
                jeq_pos = Array{Float64}(undef,tlength,1) .= 0.,
                x = Array{Float64}(undef,tlength+1,size(A,1)) .= 0.,
                y = Array{Float64}(undef,tlength,size(C,1)) .= 0.,
                Cse_Neg = Array{Float64}(undef,tlength,size(CseNegInd,1)) .= 0.,
                Cse_Pos = Array{Float64}(undef,tlength,size(CsePosInd,1)) .= 0.,
                Ce = Array{Float64}(undef,tlength,size(CeInd,1)) .= Cell.Const.ce0,
                η0 = Array{Float64}(undef,tlength,1) .= 0.,
                η_neg = Array{Float64}(undef,tlength,size(FluxNegInd,1)) .= 0.,
                ηL = Array{Float64}(undef,tlength,1) .= 0.,
                η_pos = Array{Float64}(undef,tlength,size(FluxPosInd,1)) .= 0.,
                ϕ_ẽ1 = Array{Float64}(undef,tlength,size(ϕ_ẽInd,1)) .= 0.,
                ϕ_ẽ2 = Array{Float64}(undef,tlength,size(CeInd,1)) .= 0.,
                ϕ_se_neg_0 = Array{Float64}(undef,tlength,1) .= 0., #Replace with length of ϕ_seNegInd @ zero
                jNeg = Array{Float64}(undef,tlength,size(FluxNegInd,1)) .= 0., 
                j0 = Array{Float64}(undef,tlength,1) .= 0., 
                jPos = Array{Float64}(undef,tlength,size(FluxPosInd,1)) .= 0.,
                jL = Array{Float64}(undef,tlength,1) .= 0.,
                Rtot_neg = Array{Float64}(undef,tlength,1) .= 0.,
                Rtot_pos = Array{Float64}(undef,tlength,1) .= 0., 
                Uocp_Neg = Array{Float64}(undef,tlength,1) .= 0.,
                Uocp_Pos = Array{Float64}(undef,tlength,1) .= 0.,
                Cell_V = Array{Float64}(undef,tlength,1) .= 0.,
                ϕ_e = Array{Float64}(undef,tlength,size(CeInd,1)) .= 0.,
                Cell_SOC = Array{Float64}(undef,tlength,1) .= 0,
                Iapp = Array{Float64}(undef,tlength+1,1) .= 0,
                t = t
                )


    # Defining SOC
    SOC_Neg = SOC * (Cell.Neg.θ_100-Cell.Neg.θ_0) + Cell.Neg.θ_0
    SOC_Pos = SOC * (Cell.Pos.θ_100-Cell.Pos.θ_0) + Cell.Pos.θ_0
    Results.θ_neg[1] = SOC_Neg
    Results.θ_pos[1] = SOC_Pos

    # Loop through time - compute dependent variables (voltage, flux, etc.) #
    for i in 0:(tlength-1)
        cs_neg_avg = Results.x[i+1,1] * csegain_neg + SOC_Neg * Cell.Neg.cs_max < 0. ? 0. : Results.x[i+1,1] * csegain_neg + SOC_Neg * Cell.Neg.cs_max #Zero if < 0
        cs_pos_avg = Results.x[i+1,1] * csegain_pos + SOC_Pos * Cell.Pos.cs_max < 0. ? 0. : Results.x[i+1,1] * csegain_pos + SOC_Pos * Cell.Pos.cs_max #Zero if < 0

        if cs_neg_avg > Cell.Neg.cs_max
            cs_neg_avg = Cell.Neg.cs_max
        end
        
        if cs_pos_avg > Cell.Pos.cs_max
            cs_pos_avg = Cell.Pos.cs_max
        end

        # Reaction Rates
        if Cell.Const.CellTyp == "Doyle_94"
            k_neg = Cell.Neg.k_norm
            k_pos = Cell.Pos.k_norm
            Results.jeq_neg[i+1] = k_neg*sqrt(cs_neg_avg*Cell.Const.ce0*(Cell.Neg.cs_max-cs_neg_avg))
            Results.jeq_pos[i+1] = k_pos*sqrt(cs_pos_avg*Cell.Const.ce0*(Cell.Pos.cs_max-cs_pos_avg))
        else
            k_neg = Cell.Neg.k_norm
            k_pos = Cell.Pos.k_norm
            Results.jeq_neg[i+1] = Cell.Neg.k_norm*(Cell.Const.ce0*cs_neg_avg*(Cell.Neg.cs_max-cs_neg_avg))^(1-Cell.Neg.α)
            Results.jeq_pos[i+1] = Cell.Pos.k_norm*(Cell.Const.ce0*cs_pos_avg*(Cell.Pos.cs_max-cs_pos_avg))^(1-Cell.Pos.α)
        end
        
        Results.θ_neg[i+1] = cs_neg_avg/Cell.Neg.cs_max
        Results.θ_pos[i+1] = cs_pos_avg/Cell.Pos.cs_max
        Results.Cell_SOC[i+1] = (Results.θ_neg[i+1]-Cell.Neg.θ_0)/(Cell.Neg.θ_100-Cell.Neg.θ_0)

        javg_neg = Results.Iapp[i+1]/(Cell.Neg.as*F*Cell.Neg.L*Cell.Const.CC_A)
        javg_pos = Results.Iapp[i+1]/(Cell.Pos.as*F*Cell.Pos.L*Cell.Const.CC_A)

        Arr_Factor = ((1/Cell.Const.T_ref)-(1/Tk[i+1]))/R
        κneg = Cell.Const.κf(mean(Results.Ce[i+1,CeNegInd]))*exp(Cell.Const.Ea_κ*Arr_Factor) 
        κpos = Cell.Const.κf(mean(Results.Ce[i+1,CePosInd]))*exp(Cell.Const.Ea_κ*Arr_Factor) 
        κsep = Cell.Const.κf(mean(Results.Ce[i+1,CeSepInd]))*exp(Cell.Const.Ea_κ*Arr_Factor)
        σ_neg = Cell.Neg.σ*exp(Cell.Const.Ea_κ*Arr_Factor)
        σ_pos = Cell.Pos.σ*exp(Cell.Const.Ea_κ*Arr_Factor)
        κ_eff_Neg = κneg*(Cell.Neg.ϵ_e^(Cell.Neg.κ_brug))
        κ_eff_Sep = κsep*(Cell.Sep.ϵ_e^(Cell.Sep.κ_brug))
        κ_eff_Pos = κpos*(Cell.Pos.ϵ_e^(Cell.Pos.κ_brug))
        σ_eff_Neg = σ_neg*Cell.Neg.ϵ_s^Cell.Neg.σ_brug #Effective Conductivity Neg
        σ_eff_Pos = σ_pos*Cell.Pos.ϵ_s^Cell.Pos.σ_brug #Effective Conductivity Pos
        
        # Resistances
        Results.Rtot_neg[i+1] = (Tk[i+1]*R)/(F^2*sqrt(Results.jeq_neg[i+1]^2+javg_neg^2/4))+Cell.Neg.RFilm
        Results.Rtot_pos[i+1] = (Tk[i+1]*R)/(F^2*sqrt(Results.jeq_pos[i+1]^2+javg_pos^2/4))+Cell.Pos.RFilm

        # Condensing Variable
        ν_neg = @. Cell.Neg.L*sqrt((Cell.Neg.as*(1/κ_eff_Neg+1/σ_eff_Neg))/Results.Rtot_neg[i+1])
        ν_pos = @. Cell.Pos.L*sqrt((Cell.Pos.as*(1/κ_eff_Pos+1/σ_eff_Pos))/Results.Rtot_pos[i+1])

        # Relinearise dependent on ν, σ, κ
        #D = D_Linear(Cell, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep)

        # Interpolate C & D Matrices
        C = interp(C0,SList,Results.Cell_SOC[i+1])
        D = interp(D0,SList,Results.Cell_SOC[i+1])

        # SS Output
        Results.y[i+1,:] = C*Results.x[i+1,:] + D*Results.Iapp[i+1]

        # Concentrations & Force electrode concentration maximum
        Results.Cse_Neg[i+1,:] = (SOC_Neg.*Cell.Neg.cs_max .+ Results.y[i+1,CseNegInd]) > ones(size(Results.Cse_Neg,2))*Cell.Neg.cs_max ? ones(size(Results.Cse_Neg,2))*Cell.Neg.cs_max : (SOC_Neg.*Cell.Neg.cs_max .+ Results.y[i+1,CseNegInd])
        Results.Cse_Pos[i+1,:] = (SOC_Pos.*Cell.Pos.cs_max .+ Results.y[i+1,CsePosInd]) > ones(size(Results.Cse_Pos,2))*Cell.Pos.cs_max ? ones(size(Results.Cse_Pos,2))*Cell.Pos.cs_max : (SOC_Pos.*Cell.Pos.cs_max .+ Results.y[i+1,CsePosInd])
        Results.Ce[i+1,:] = @. Cell.Const.ce0 + Results.y[i+1,CeInd]
    
        # Potentials
        Results.Uocp_Neg[i+1] = Cell.Const.Uocp("Neg",Results.Cse_Neg[i+1,1]/Cell.Neg.cs_max)
        Results.Uocp_Pos[i+1] = Cell.Const.Uocp("Pos",Results.Cse_Pos[i+1,1]/Cell.Pos.cs_max)
        Results.ϕ_se_neg_0[i+1] = Results.y[i+1, ϕ_seNegInd[1]] + Results.Uocp_Neg[i+1] #Location 0
        Results.ϕ_ẽ1[i+1,:] = Results.y[i+1,ϕ_ẽInd]
        Results.ϕ_ẽ2[i+1,:] = @. ((Tk[i+1]*2*R*(1-Cell.Const.tpf(Results.Ce[i+1,:])))/F)*(log(Results.Ce[i+1,:]/Results.Ce[i+1,1]))
        Results.ϕ_e[i+1,:] = @. [0; Results.ϕ_ẽ1[i+1,:]]+Results.ϕ_ẽ2[i+1,:]-Results.ϕ_se_neg_0[i+1]

        # Flux
        Results.jNeg[i+1,:] = Results.y[i+1,FluxNegInd]
        Results.j0[i+1] = Results.y[i+1,FluxNegInd[1]]
        Results.jPos[i+1,:] = Results.y[i+1,FluxPosInd]
        Results.jL[i+1] = Results.y[i+1,FluxPosInd[1]]

            # Neg
            j0_neg = findmax([ones(size(Results.Cse_Neg,2))*eps() (k_neg.*(Results.Cse_Neg[i+1,:].^Cell.Neg.α).*(Results.Ce[i+1,1].^(1-Cell.Neg.α))).*(Cell.Neg.cs_max.-Results.Cse_Neg[i+1,:]).^(1-Cell.Neg.α)], dims=2)[1]
            Results.η0[i+1] = Tk[i+1]*2*R/F*asinh(Results.j0[i+1]/(2*j0_neg[1]))
            Results.η_neg[i+1,:] = @. (Tk[i+1]*2*R)/F*asinh((Results.jNeg[i+1,:])/(2*j0_neg))

            # Pos
            j0_pos = findmax([ones(size(Results.Cse_Pos,2))*eps() (k_pos.*(Results.Cse_Pos[i+1,:].^Cell.Pos.α).*(Results.Ce[i+1,1].^(1-Cell.Pos.α))).*(Cell.Pos.cs_max.-Results.Cse_Pos[i+1,:]).^(1-Cell.Pos.α)], dims=2)[1]
            Results.ηL[i+1] = (Tk[i+1]*2*R)/F*asinh(Results.jL[i+1]/(2*j0_pos[1]))
            Results.η_pos[i+1,:] = @. (Tk[i+1]*2*R)/F*asinh(Results.jPos[i+1,:]/(2*j0_pos))
    

        # Cell Voltage
        Results.Cell_V[i+1] = @. (Results.Uocp_Pos[i+1]-Results.Uocp_Neg[i+1]) + (Results.ηL[i+1]-Results.η0[i+1]) + (Results.ϕ_ẽ1[i+1,end]+Results.ϕ_ẽ2[i+1,end]) + (Cell.Pos.RFilm*Results.jL[i+1]-Cell.Neg.RFilm*Results.j0[i+1])*F 

        # ϕ_s
        ϕ_s_neg = Results.y[i+1,ϕ_sNegInd]
        ϕ_s_pos = @. Results.y[i+1,ϕ_sPosInd] + Results.Cell_V[i+1]

        # Interpolate A Matrix
        A = interp(A0,SList,Results.Cell_SOC[i+1])
        B = interp(B0,SList,Results.Cell_SOC[i+1])

        # Update States
        Results.x[i+2,:] = A*Results.x[i+1,:] + B*Results.Iapp[i+1]


        if Def == "Power"
            Results.Iapp[i+2] = Input[i+1,1]/Cell_V[i+1]
        else
            Results.Iapp[i+2] = Input[i+1,1]
        end

    end

return Results
end