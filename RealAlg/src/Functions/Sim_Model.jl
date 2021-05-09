function Sim_Model(CellData,Iapp,Tk,A,B,C,D)
    """ 
    Simulation of generated reduced-order models
    # Add License
    # Add Ins and Outs
        # Locations to be computed
        # Sampling Frequency
    """

    #Determine time span and allocate arrays
    tlength = size(Iapp,1)
    #Determine Inital conditions

    #Defining SOC
    θ_neg[1] = CellData.Const.SOC * (CellData.Neg.θ_100-CellData.Neg.θ_0) + CellData.Neg.θ_0
    θ_pos[1] = CellData.Const.SOC * (CellData.Pos.θ_100-CellData.Pos.θ_0) + CellData.Pos.θ_0
    SOC = CellData.Concentration.SOC




    #Loop through time
        #Compute dependent variables (voltage, flux, etc.)


    for i in 2:tlength
        
        cs_neg_avg = x[k+1,end]*csegain_neg+θ_neg*CellData.Neg.cs_max
        θ_neg = cs_neg_avg/CellData.Neg.cs_max
        cs_pos_avg = x[k+1,end]*csegain_pos+θ_pos*CellData.Pos.cs_max
        θ_pos = cs_pos_avg/CellData.Pos.cs_max
        Cell_SOC = (θ_neg-CellData.Neg.θ_0)/(CellData.Neg.θ_100-CellData.Neg.θ_0)

        jeq_neg = CellData.Neg.k*sqrt(cs_neg_avg*(CellData.Const.Ce0*(CellData.Neg.cs_max-cs_neg_avg)))
        jeq_pos = CellData.Pos.k*sqrt(cs_pos_avg*(CellData.Const.Ce0*(CellData.Pos.cs_max-cs_pos_avg)))
        javg_neg = Iapp[i]/(CellData.Neg.ϵ_s*(3*F*CellData.Neg.L*CellData.Neg.CC_A)/CellData.Neg.Rs)
        javg_neg = Iapp[i]/(CellData.Neg.ϵ_s*(3*F*CellData.Pos.L*CellData.Pos.CC_A)/CellData.Pos.Rs)


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
        k_norm_neg = CellData.Neg.k*CellData.Neg.cs_max/CellData.Const.Ce0^(1-CellData.Neg.α)
        k_norm_pos = CellData.Pos.k*CellData.Pos.cs_max/CellData.Const.Ce0^(1-CellData.Pos.α)

        #Resistances
        Rtot_neg = (Tk[i]*R)/(F^2*sqrt(jeq_neg^2+(javg_neg^2/4)))+CellData.Neg.RFilm
        Rtot_pos = (Tk[i]*R)/(F^2*sqrt(jeq_pos^2+(javg_pos^2/4)))+CellData.Pos.RFilm

        #Condensing Variable
        ν_neg = CellData.Neg.L*sqrt((3*(CellData.Neg.ϵ_s/CellData.Neg.Rs)*(1/κ_eff_Neg+1/σ_eff_Neg))/Rtot_neg)
        ν_pos = CellData.Pos.L*sqrt((3*(CellData.Pos.ϵ_s/CellData.Pos.Rs)*(1/κ_eff_Pos+1/σ_eff_Pos))/Rtot_pos)


        #Relinearise dependent on ν, σ, κ
        #Call from CellData? List of functions composed from ROM creation?
        D = D_fun(CellData, ν_neg, ν_pos, σ_eff_Neg, σ_eff_Pos, κ_eff_Neg, κ_eff_Sep, κ_eff_Pos) #Calling D linearisation functions


        y[i,:] = C*x[i,:]+D*Iapp[i]

        
        #Concentrations
        Cse_Neg_0 = 
        Cse_Neg_1 = 
        Cse_Pos_0 =
        Cse_Pos_1 = 
        Ce_0 = 
        Ce_1 = 
        Ce_2 = 
        Ce_3 = 

        #Potentials
        Phi_se_0 = 
        Uocp_Neg = 
        Uocp_Pos = 
        Phi_e_tilde1 = 
        Phi_e_tilde2
        Phi_e = 

        #Flux
        j0_CC_neg = 
        j0_neg =  
        η0 = 
        η_neg = 
        j0_CC_pos = 
        j0_pos = 
        ηL = 
        η_pos = 

        #Cell Voltage
        Cell_V = 
        η_overP_neg = 

        #
        x[i+1,:] = A*x[i,:]'+B*Iapp[i]

    end
    #Join data into struct and export

end
