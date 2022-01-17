function Sim_Model(CellData,Iapp,Tk,SOC,A0,B0,C0,D0)
    """ 
    Simulation of generated reduced-order models
    # Add License
    # Add Ins and Outs
        # Locations to be computed
        # Sampling Frequency
    """
    #Determine time span and allocate arrays
    tlength = size(Iapp,1)

    #CellV_= Array{Float64}(undef,tlength,0)
    Ce_= Array{Float64}(undef,tlength,0)
    j0_ = Array{Float64}(undef,tlength,0)
    RtotNeg_ = Array{Float64}(undef,tlength,0)
    RtotPos_ = Array{Float64}(undef,tlength,0)
    ϕ_ẽ1_ = Array{Float64}(undef,tlength,0)
    ϕ_ẽ2_ = Array{Float64}(undef,tlength,0)
    Uocp_Neg_ = Array{Float64}(undef,tlength,0)
    Uocp_Pos_ = Array{Float64}(undef,tlength,0)
    η0_ = Array{Float64}(undef,tlength,0)
    η_neg_ = Array{Float64}(undef,tlength,0)
    ηL_ = Array{Float64}(undef,tlength,0)
    η_pos_ = Array{Float64}(undef,tlength,0)
    ϕ_e_ = Array{Float64}(undef,tlength,0)
    jNeg_ = Array{Float64}(undef,tlength,0)
    jPos_ = Array{Float64}(undef,tlength,0)
    CellV_ = Cse_Neg_ = Cse_Pos_ = tuple()
    #Cse_Pos_ = Array{Float64}(undef,tlength,0)

    #Selecting SS Models
    for γ in 1:1:tuple_len(A0)
        A = A0[γ]
        B = B0[γ]
        C = C0[γ]
        D = D0[γ]
        #@show CellData.Const.SOC = SOC[γ]
        CellData.Const.SOC = SOC
        #Capturing Indices
        tfstr = Array{String}(undef,0,1)
        for i in 1:size(CellData.Transfer.tfs[:,1],1)
            t1 = t2 = Array{String}(undef,0,1)
            for j in 1:size(CellData.Transfer.tfs[i,3],1)
                t1 ="$(CellData.Transfer.tfs[i,1])_$(CellData.Transfer.tfs[i,2])"
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

        CeNeg = CellData.Transfer.tfs[1,3].<=CellData.Neg.L
        CeNegOffset = CellData.Transfer.tfs[1,3].<CellData.Neg.L
        CeSep = CellData.Transfer.tfs[1,3].<=CellData.Neg.L+CellData.Sep.L
        CeSepOffset = CellData.Transfer.tfs[1,3].<CellData.Neg.L+CellData.Sep.L
        CePos = CellData.Transfer.tfs[1,3].<=CellData.Const.Ltot
        CeNegInd = findall(CeNeg .== 1)
        CeSepInd = findall(CeSep.-CeNegOffset .== 1)
        CePosInd = findall(CePos.-CeSepOffset .== 1)

        csegain_neg = C[CseNegInd[1][1],end]
        csegain_pos = C[CsePosInd[1][1],end]

        #Memory Allocation
        θ_neg = Array{Float64}(undef,tlength,1) .= 0.
        θ_pos = Array{Float64}(undef,tlength,1) .= 0.
        x = Array{Float64}(undef,tlength,size(A,1)) .= 0.
        y = Array{Float64}(undef,tlength,size(C,1)) .= 0.
        Cse_Neg = Array{Float64}(undef,tlength,size(CseNegInd,1)) .= 0.
        Cse_Pos = Array{Float64}(undef,tlength,size(CsePosInd,1)) .= 0.
        Ce = Array{Float64}(undef,tlength,size(CeInd,1)) .= CellData.Const.ce0
        η0 = Array{Float64}(undef,tlength,1) .= 0.
        η_neg = Array{Float64}(undef,tlength,size(FluxNegInd,1)) .= 0.
        ηL = Array{Float64}(undef,tlength,1) .= 0.
        η_pos = Array{Float64}(undef,tlength,size(FluxPosInd,1)) .= 0.
        ϕ_ẽ1 = Array{Float64}(undef,tlength,size(ϕ_ẽInd,1)) .= 0.
        ϕ_ẽ2 = Array{Float64}(undef,tlength,size(CeInd,1)) .= 0.
        ϕ_se_neg_0 = Array{Float64}(undef,tlength,1) .= 0. #Replace with length of ϕ_seNegInd @ zero
        jNeg = Array{Float64}(undef,tlength,size(FluxNegInd,1)) .= 0. 
        j0 = Array{Float64}(undef,tlength,1) .= 0. 
        jPos = Array{Float64}(undef,tlength,size(FluxPosInd,1)) .= 0.
        jL = Array{Float64}(undef,tlength,1) .= 0. 
        j0_CC_neg = Array{Float64}(undef,tlength,1) .= 0.
        Rtot_neg = Array{Float64}(undef,tlength,1) .= 0.
        Rtot_pos = Array{Float64}(undef,tlength,1) .= 0.
        j0_CC_pos = Array{Float64}(undef,tlength,1) .= 0.  
        Uocp_Neg = Array{Float64}(undef,tlength,1) .= 0.
        Uocp_Pos = Array{Float64}(undef,tlength,1) .= 0.
        Cell_V = Array{Float64}(undef,tlength,1) .= 0.
        ϕ_e = Array{Float64}(undef,tlength,size(CeInd,1)) .= 0.


        #Defining SOC
        SOC_Neg = CellData.Const.SOC * (CellData.Neg.θ_100-CellData.Neg.θ_0) + CellData.Neg.θ_0
        SOC_Pos = CellData.Const.SOC * (CellData.Pos.θ_100-CellData.Pos.θ_0) + CellData.Pos.θ_0
        θ_neg[1] = SOC_Neg
        θ_pos[1] = SOC_Pos

        #Loop through time
        #Compute dependent variables (voltage, flux, etc.)
        for i in 1:tlength-1
            cs_neg_avg = x[i,end] * csegain_neg + SOC_Neg * CellData.Neg.cs_max < 0. ? 0. : x[i,end] * csegain_neg + SOC_Neg * CellData.Neg.cs_max #Zero if < 0
            cs_pos_avg = x[i,end] * csegain_pos + SOC_Pos * CellData.Pos.cs_max < 0. ? 0. : x[i,end] * csegain_pos + SOC_Pos * CellData.Pos.cs_max #Zero if < 0

            if cs_neg_avg > CellData.Neg.cs_max
                cs_neg_avg = CellData.Neg.cs_max
            end
            
            if cs_pos_avg > CellData.Pos.cs_max
                cs_pos_avg = CellData.Pos.cs_max
            end

            #Reaction Rates
            if CellData.Const.CellTyp == "Doyle_94"
                k_neg = CellData.Neg.k_norm/CellData.Neg.cs_max/CellData.Const.ce0^(1-CellData.Neg.α)
                k_pos = CellData.Pos.k_norm/CellData.Pos.cs_max/CellData.Const.ce0^(1-CellData.Pos.α)
                jeq_neg = k_neg*sqrt(cs_neg_avg*CellData.Const.ce0*(CellData.Neg.cs_max-cs_neg_avg))
                jeq_pos = k_pos*sqrt(cs_pos_avg*CellData.Const.ce0*(CellData.Pos.cs_max-cs_pos_avg))
            else
                jeq_neg = CellData.Neg.k_norm*CellData.Neg.cs_max*(CellData.Const.ce0*(cs_neg_avg/CellData.Neg.cs_max*(1-cs_neg_avg/CellData.Neg.cs_max)))^(1-CellData.Neg.α)
                jeq_pos = CellData.Pos.k_norm*CellData.Pos.cs_max*(CellData.Const.ce0*(cs_pos_avg/CellData.Pos.cs_max*(1-cs_pos_avg/CellData.Pos.cs_max)))^(1-CellData.Pos.α)
            end
            
            θ_neg[i+1] = cs_neg_avg/CellData.Neg.cs_max
            θ_pos[i+1] = cs_pos_avg/CellData.Pos.cs_max
            Cell_SOC = (θ_neg[i]-CellData.Neg.θ_0)/(CellData.Neg.θ_100-CellData.Neg.θ_0)


            javg_neg = Iapp[i]/(CellData.Neg.as*F*CellData.Neg.L*CellData.Const.CC_A)
            javg_pos = Iapp[i]/(CellData.Pos.as*F*CellData.Pos.L*CellData.Const.CC_A)

            Arr_Factor = ((1/CellData.Const.T_ref)-(1/Tk[i]))/R
            κneg = CellData.Const.κf(mean(Ce[i,CeNegInd]))*exp(CellData.Const.Ea_κ*Arr_Factor) 
            κpos = CellData.Const.κf(mean(Ce[i,CePosInd]))*exp(CellData.Const.Ea_κ*Arr_Factor) 
            κsep = CellData.Const.κf(mean(Ce[i,CeSepInd]))*exp(CellData.Const.Ea_κ*Arr_Factor)
            σ_neg = CellData.Neg.σ*exp(CellData.Const.Ea_κ*Arr_Factor)
            σ_pos = CellData.Pos.σ*exp(CellData.Const.Ea_κ*Arr_Factor)
            κ_eff_Neg = κneg*(CellData.Neg.ϵ_e^(CellData.Neg.κ_brug))
            κ_eff_Sep = κsep*(CellData.Sep.ϵ_e^(CellData.Sep.κ_brug))
            κ_eff_Pos = κpos*(CellData.Pos.ϵ_e^(CellData.Pos.κ_brug))
            σ_eff_Neg = σ_neg*CellData.Neg.ϵ_s^CellData.Neg.σ_brug #Effective Conductivity Neg
            σ_eff_Pos = σ_pos*CellData.Pos.ϵ_s^CellData.Pos.σ_brug #Effective Conductivity Pos
            
            #Resistances
            Rtot_neg[i] = (Tk[i]*R)/(F^2*mean([abs(jeq_neg) abs(javg_neg)]))+CellData.Neg.RFilm
            Rtot_pos[i] = (Tk[i]*R)/(F^2*mean([abs(jeq_pos) abs(javg_pos)]))+CellData.Pos.RFilm

            #Condensing Variable
            ν_neg = @. CellData.Neg.L*sqrt((CellData.Neg.as*(1/κ_eff_Neg+1/σ_eff_Neg))/Rtot_neg[i])
            ν_pos = @. CellData.Pos.L*sqrt((CellData.Pos.as*(1/κ_eff_Pos+1/σ_eff_Pos))/Rtot_pos[i])

            #Relinearise dependent on ν, σ, κ
            #Call from CellData? List of functions composed from ROM creation?
            D = D_Linear(CellData, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep)

            #SS Output
            y[i,:] = C*x[i,:] + D*Iapp[i]

            #Concentrations & Force electrode concentration maximum
            Cse_Neg[i,:] = (SOC_Neg.*CellData.Neg.cs_max .+ y[i,CseNegInd]) > ones(size(Cse_Neg,2))*CellData.Neg.cs_max ? ones(size(Cse_Neg,2))*CellData.Neg.cs_max : (SOC_Neg.*CellData.Neg.cs_max .+ y[i,CseNegInd])
            Cse_Pos[i,:] = (SOC_Pos.*CellData.Pos.cs_max .+ y[i,CsePosInd]) > ones(size(Cse_Pos,2))*CellData.Pos.cs_max ? ones(size(Cse_Pos,2))*CellData.Pos.cs_max : (SOC_Pos.*CellData.Pos.cs_max .+ y[i,CsePosInd])
            Ce[i,:] = @. CellData.Const.ce0 + y[i,CeInd]
        
            #Potentials
            Uocp_Neg[i] = CellData.Const.Uocp("Neg",Cse_Neg[i,1]/CellData.Neg.cs_max)
            Uocp_Pos[i] = CellData.Const.Uocp("Pos",Cse_Pos[i,1]/CellData.Pos.cs_max)
            ϕ_se_neg_0[i] = y[i, ϕ_seNegInd[1]] + Uocp_Neg[i] #Location 0
            ϕ_ẽ1[i,:] = y[i,ϕ_ẽInd]
            ϕ_ẽ2[i,:] = @. ((Tk[i]*2*R*(1-CellData.Const.tpf(Ce[i,:])))/F)*(log(Ce[i,:]/Ce[i,1]))
            ϕ_e[i,:] = @. [0; ϕ_ẽ1[i,:]]+ϕ_ẽ2[i,:]-ϕ_se_neg_0[i]

            #Flux
            jNeg[i,:] = y[i,FluxNegInd]
            j0[i] = y[i,FluxNegInd[1]]
            jPos[i,:] = y[i,FluxPosInd]
            jL[i] = y[i,FluxPosInd[1]]


            j0_CC_neg[i] = findmax([eps(); (k_neg*(CellData.Neg.cs_max-Cse_Neg[i,1])^(1-CellData.Neg.α))*((Cse_Neg[i,1]^CellData.Neg.α)*(Ce[i,1]^(1-CellData.Neg.α)))])[1]
            j0_neg = findmax([ones(size(Cse_Neg,2))*eps() (k_neg.*(Cse_Neg[i,:].^CellData.Neg.α).*(Ce[i,1].^(1-CellData.Neg.α))).*(CellData.Neg.cs_max.-Cse_Neg[i,:]).^(1-CellData.Neg.α)], dims=2)[1]
            η0[i] = Tk[i]*2*R/F*asinh(j0[i]/(2*j0_CC_neg[i]))
            η_neg[i,:] = @. (Tk[i]*2*R)/F*asinh((jNeg[i,:])/(2*j0_neg))

            j0_CC_pos[i] = findmax([eps(); (k_pos*(CellData.Pos.cs_max-Cse_Pos[i,1])^(1-CellData.Pos.α))*((Cse_Pos[i,1]^CellData.Pos.α)*(Ce[i,end]^(1-CellData.Pos.α)))])[1] 
            j0_pos = findmax([ones(size(Cse_Pos,2))*eps() (k_pos.*(Cse_Pos[i,:].^CellData.Pos.α).*(Ce[i,1].^(1-CellData.Pos.α))).*(CellData.Pos.cs_max.-Cse_Pos[i,:]).^(1-CellData.Pos.α)], dims=2)[1]
            ηL[i] = (Tk[i]*2*R)/F*asinh(jL[i]/(2*j0_CC_pos[i]))
            η_pos[i,:] = @. (Tk[i]*2*R)/F*asinh(jPos[i,:]/(2*j0_pos))
        

            #Cell Voltage
            Cell_V[i] = @. (Uocp_Pos[i]-Uocp_Neg[i]) + (ηL[i]-η0[i]) + (ϕ_ẽ1[i,end]+ϕ_ẽ2[i,end]) + (CellData.Pos.RFilm*jL[i]-CellData.Neg.RFilm*j0[i])*F 

            #ϕ_s
            ϕ_s_neg = y[i,ϕ_sNegInd]
            ϕ_s_pos = @. y[i,ϕ_sPosInd] + Cell_V[i]

            #Update States
            x[i+1,:] = A*x[i,:] + B*Iapp[i]
        end

        CellV_ = flatten_(CellV_, Cell_V)
        Ce_ = [Ce_ Ce]
        jNeg_ = jNeg
        jPos_ = jPos
        RtotNeg_ = [RtotNeg_ Rtot_neg]
        RtotPos_ = [RtotPos_ Rtot_pos]
        η0_ = [η0_ η0]
        η_neg_ = η_neg
        ηL_ = [ηL_ ηL]
        η_pos_ = η_pos
        ϕ_ẽ1_ = ϕ_ẽ1
        ϕ_ẽ2_ = ϕ_ẽ2
        ϕ_e_ = ϕ_e
        Uocp_Neg_ = [Uocp_Neg_ Uocp_Neg]
        Uocp_Pos_ = [Uocp_Pos_ Uocp_Pos]
        Cse_Neg_ = flatten_(Cse_Neg_, Cse_Neg)
        Cse_Pos_ = flatten_(Cse_Pos_, Cse_Pos)
    end 
    return CellV_, Ce_, jNeg_, jPos_, RtotNeg_, RtotPos_, η0_, ηL_, η_neg_, η_pos_, ϕ_ẽ1_, ϕ_ẽ2_, Uocp_Neg_, Uocp_Pos_, ϕ_e_, Cse_Neg_, Cse_Pos_
end