function Sim_Model(CellData,Dtt,Iapp,Tk,A,B,C,D)
    """ 
    Simulation of generated reduced-order models
    # Add License
    # Add Ins and Outs
        # Locations to be computed
        # Sampling Frequency
    """
    #Slice SS Matrices
    A = A[1]
    B = B[1]
    C = C[1]
    D = D[1]

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
    PhisNegInd = findall(isequal("Phi_s_Neg"), tfstr)
    PhisPosInd = findall(isequal("Phi_s_Pos"), tfstr)
    ϕ_sePosInd = findall(isequal("Phi_se_Pos"), tfstr)
    ϕ_seNegInd = findall(isequal("Phi_se_Neg"), tfstr)
    FluxNegInd = findall(isequal("Flux_Neg"), tfstr)
    FluxPosInd = findall(isequal("Flux_Pos"), tfstr)

    csegain_neg = C[CseNegInd[1],1]
    csegain_pos = C[CsePosInd[1],1]

    #Determine time span and allocate arrays
    tlength = size(Iapp,1)

    #Memory Allocation
    θ_neg = θ_pos = Array{Float64}(undef,tlength,1)
    x = Array{Float64}(undef,1,size(A,1)) .= 0.
    y = Array{Float64}(undef,1,size(C,1)) .= 0.
    Cse_Neg = Array{Float64}(undef,1,size(CseNegInd,1)) .= 0.
    Cse_Pos = Array{Float64}(undef,1,size(CsePosInd,1)) .= 0.
    Ce = Array{Float64}(undef,1,size(CeInd,1)) .= 0.
    η0 = Array{Float64}(undef,1,size(CeInd,1)) .= 0.
    η_neg = Array{Float64}(undef,1,size(FluxNegInd,1)) .= 0.
    ηL = Array{Float64}(undef,1,size(CeInd,1)) .= 0.
    η_pos = Array{Float64}(undef,1,size(FluxPosInd,1)) .= 0.


    #Defining SOC
    θ_neg[1] = CellData.Const.SOC * (CellData.Neg.θ_100-CellData.Neg.θ_0) + CellData.Neg.θ_0
    θ_pos[1] = CellData.Const.SOC * (CellData.Pos.θ_100-CellData.Pos.θ_0) + CellData.Pos.θ_0
    SOC = CellData.Const.SOC

    #Loop through time
        #Compute dependent variables (voltage, flux, etc.)
    for i in 1:tlength
        cs_neg_avg = x[i,end]*csegain_neg+θ_neg[i]*CellData.Neg.cs_max < 0. ? 0. : x[i,end]*csegain_neg+θ_neg[i]*CellData.Neg.cs_max #Zero if < 0
        cs_pos_avg = x[i,end]*csegain_pos+θ_pos[i]*CellData.Pos.cs_max < 0. ? 0. : x[i,end]*csegain_pos+θ_pos[i]*CellData.Pos.cs_max #Zero if < 0

        println("cs_neg_avg:", cs_neg_avg)
        θ_neg[i] = cs_neg_avg/CellData.Neg.cs_max

        θ_pos[i] = cs_pos_avg/CellData.Pos.cs_max
        Cell_SOC = (θ_neg[i]-CellData.Neg.θ_0)/(CellData.Neg.θ_100-CellData.Neg.θ_0)

        jeq_neg = CellData.Neg.k_norm*sqrt(cs_neg_avg*(CellData.Const.ce0*(CellData.Neg.cs_max-cs_neg_avg)))
        jeq_pos = CellData.Pos.k_norm*sqrt(cs_pos_avg*(CellData.Const.ce0*(CellData.Pos.cs_max-cs_pos_avg)))
        javg_neg = Iapp[i]/(CellData.Neg.ϵ_s*(3*F*CellData.Neg.L*CellData.Const.CC_A)/CellData.Neg.Rs)
        javg_pos = Iapp[i]/(CellData.Pos.ϵ_s*(3*F*CellData.Pos.L*CellData.Const.CC_A)/CellData.Pos.Rs)


        Arr_Factor = ((1/CellData.Const.T_ref)-(1/Tk[i]))/R
        κ = CellData.Const.κf(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
        σ_neg = CellData.Neg.σ*exp(CellData.Const.Ea_κ*Arr_Factor)
        σ_pos = CellData.Pos.σ*exp(CellData.Const.Ea_κ*Arr_Factor)
        κ_eff_Neg = κ*(CellData.Neg.ϵ_e^(CellData.Neg.κ_brug))
        κ_eff_Sep = κ*(CellData.Sep.ϵ_e^(CellData.Sep.κ_brug))
        κ_eff_Pos = κ*(CellData.Pos.ϵ_e^(CellData.Pos.κ_brug))


        σ_eff_Neg = σ_neg*CellData.Neg.ϵ_s^CellData.Neg.σ_brug #Effective Conductivity Neg
        σ_eff_Pos = σ_pos*CellData.Pos.ϵ_s^CellData.Pos.σ_brug #Effective Conductivity Pos
        
        #Reaction Rates
        k_neg = CellData.Neg.k_norm*CellData.Neg.cs_max/CellData.Const.ce0^(1-CellData.Neg.α)
        k_pos = CellData.Pos.k_norm*CellData.Pos.cs_max/CellData.Const.ce0^(1-CellData.Pos.α)

        #Resistances
        Rtot_neg = (Tk[i]*R)/(F^2*sqrt(jeq_neg^2+(javg_neg^2/4)))+CellData.Neg.RFilm
        Rtot_pos = (Tk[i]*R)/(F^2*sqrt(jeq_pos^2+(javg_pos^2/4)))+CellData.Pos.RFilm

        #Condensing Variable
        ν_neg = CellData.Neg.L*sqrt((3*(CellData.Neg.ϵ_s/CellData.Neg.Rs)*(1/κ_eff_Neg+1/σ_eff_Neg))/Rtot_neg)
        ν_pos = CellData.Pos.L*sqrt((3*(CellData.Pos.ϵ_s/CellData.Pos.Rs)*(1/κ_eff_Pos+1/σ_eff_Pos))/Rtot_pos)


        #Relinearise dependent on ν, σ, κ
        #Call from CellData? List of functions composed from ROM creation?
        #D = D_fun(CellData, ν_neg, ν_pos, σ_eff_Neg, σ_eff_Pos, κ_eff_Neg, κ_eff_Sep, κ_eff_Pos) #Calling D linearisation functions
        Dtemp = Array{Float64}(undef,0,1)
        D = Array{Float64}(undef,0,1)
        for j in 1:length(Dtt)
            Dtemp = eval(Meta.parse(Dtt[j]))
            D = [D; Dtemp]
        end

        #SS Output
        y[i,:] = C*x[i,:]+D*Iapp[i]
        println("y:",size(y))
        display("text/plain", y[i,:])

        #Concentrations
        Cse_Neg[i,:] = @. SOC*CellData.Neg.cs_max + y[i,CseNegInd]  
        Cse_Pos[i,:] = @. SOC*CellData.Pos.cs_max + y[i,CsePosInd] 
        Ce[i,:] = @. CellData.Const.ce0 + y[i,CeInd]

        #Potentials
        Uocp_Neg = CellData.Const.Uocp("Neg",θ_neg)
        Uocp_Pos = CellData.Const.Uocp("Pos",θ_pos)
        ϕ_se_neg_0 = y[i, ϕ_seNegInd[1]] #Location 0
        ϕ_ẽ1 = [0; y[i,ϕ_ẽInd]]
        ϕ_ẽ2 = @. ((Tk[i]*2*R)/F)*(log(Ce/Ce[1])*(1-CellData.Const.t_plus))
        ϕ_e = @. ϕ_ẽ1+ϕ_ẽ2-ϕ_se_neg_0

        #Flux
        j0_CC_neg = findmax(((CellData.Neg.cs_max+Cse_Neg[i,1])^(1-CellData.Neg.α))*((Cse_Neg[i,1]^CellData.Neg.α)*(Ce[1]^(1-CellData.Neg.α)))*k_neg)[1]
        j0_neg = @. findmax(((Cse_Neg^CellData.Neg.α)*(Ce[i,1:2]^(1-CellData.Neg.α)))*(CellData.Neg.cs_max-Cse_Neg[i,:])^(1-CellData.Neg.α)*k_neg)[1]
        η0[i,:] = @. asinh((y[i,FluxNegInd[1]]/(2*j0_CC_neg)))*(Tk[i]*2*R/F)
        η_neg[i,:] = @. (Tk[i]*2*R/F)/asinh(y[i,FluxNegInd]/(2*j0_neg[i,:]))

        j0_CC_pos = findmax(((CellData.Pos.cs_max+Cse_Pos[1])^(1-CellData.Pos.α))*((Cse_Pos[1]^CellData.Pos.α)*(Ce[1]^(1-CellData.Pos.α)))*k_pos)[1] 
        j0_pos = @. findmax(((Cse_Pos^CellData.Pos.α)*(Ce[1:2]^(1-CellData.Pos.α)))*(CellData.Pos.cs_max-Cse_Pos)^(1-CellData.Pos.α)*k_pos)[1]
        ηL[i,:] = @. asinh((y[i,FluxPosInd[1]]/(2*j0_CC_pos)))*(Tk[i]*2*R/F)
        η_po[i,:] = @. (Tk[i]*2*R/F)/asinh((y[i,FluxPosInd]/(2*j0_pos[i,:])))

        #Cell Voltage
        Cell_V[i] = (Uocp_Pos-Uocp_Neg) + (ηL-η0) + ϕ_ẽ1[end] + ϕ_ẽ2[end]+ (CellData.Pos.RFilm*j0_pos-CellData.Neg.RFilm*j0_neg)*F

        #ϕ_s
        ϕ_s_neg = y[i,ϕ_sNegInd]
        ϕ_s_pos = y[i,ϕ_sPosInd] + Cell_V[i]


        #Update States
        x[i+1,:] = A*x[i,:]' + B*Iapp[i]

    end
    return Cell_V
end
