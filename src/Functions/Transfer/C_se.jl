@inline function C_se(CellData,s,z,Def)
    """ 
    Concentration Solid-Electrolyte Transfer Function
    # Add Ins and Outs
        # Cell Data 
        # Frequency Vector 
        # Discretisation Locations
        # Electrode Definition
    """


 if Def == "Pos"   
    Electrode = CellData.Pos #Electrode Length
 else
    Electrode = CellData.Neg #Electrode Length
 end


κ_eff = CellData.Const.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = CellData.Const.SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0 

#Beta's
β = @. Electrode.Rs*sqrt(s/Electrode.Ds)

#Prepare for j0
cs0 = Electrode.cs_max * θ

#Current Flux Density
if CellData.Const.CellTyp == "Doyle_94"
    κ = Electrode.k_norm/Electrode.cs_max/CellData.Const.ce0^(1-Electrode.α)
    j0 = κ*(CellData.Const.ce0*(Electrode.cs_max-cs0))^(1-Electrode.α)*cs0^Electrode.α
else
    j0 = Electrode.k_norm*(CellData.Const.ce0*(cs0/Electrode.cs_max*(1-cs0/Electrode.cs_max)))^(1-Electrode.α)
end

#Resistances
Rtot = R*CellData.Const.T/(j0*F^2) + Electrode.RFilm

#∂Uocp_Def
∂Uocp_elc = CellData.Const.∂Uocp(Def,θ)/Electrode.cs_max

cse_res = -3/(Electrode.as*F*Electrode.L*CellData.Const.CC_A*Electrode.Rs) #Residual Variable for Pole Removal - eq. 4.45
ν = @. Electrode.L*sqrt((Electrode.as/σ_eff+Electrode.as/κ_eff)/(Rtot+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(β)/(tanh(β)-β)))) #Condensing Variable - eq. 4.13

cse_tf = @. ν*Electrode.Rs*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))*tanh(β)/(Electrode.as*F*Electrode.L*CellData.Const.CC_A*Electrode.Ds*(κ_eff+σ_eff)*sinh(ν)*(tanh(β)-β)) #Transfer Function - eq. 4.17
cse_tf = @. cse_tf-cse_res/s  #Pole removal - eq. 4.43
zero_tf = @. (5*Electrode.as*Electrode.Ds*F*Electrode.L^2*(κ_eff*(2-6*z+3*z^2)+(3*z^2-1)*σ_eff)-6*∂Uocp_elc*Electrode.Rs*κ_eff*σ_eff)/(30*CellData.Const.CC_A*Electrode.as*Electrode.Ds*∂Uocp_elc*F*Electrode.L*κ_eff*σ_eff) #For s = 0 / Wolfram Alpha
cse_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]
D = zeros(length(z))
D_term = "zeros(length($z))"

if Def == "Pos"
    cse_tf = -cse_tf
    cse_res = -cse_res

    if Debug == 1
        println("D:C_se:Pos",D)
        println("z:C_se:Pos",z)
        println("zero_tf:C_se:Pos",zero_tf[:,1])
        println("κ_eff:C_se:Pos",κ_eff)
        println("σ_eff:C_se:Pos",σ_eff)
        #println("ν_∞:C_se:Pos",ν_∞)
        println("j0:C_se:Pos",j0)
        println("as:C_se:Pos",Electrode.as)
        println("Rtot:C_se:Pos",Rtot)
        println("σ_eff:C_se:Pos",σ_eff)
        println("∂Uocp_elc:C_se:Pos",∂Uocp_elc)
        println("Electrode.L:C_se:Pos",Electrode.L)
    end
end
cse_res = cse_res*ones(length(z))
return cse_tf, D, cse_res, D_term

end
