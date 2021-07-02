@inline function C_se(CellData,s,z,Def)
    """ 
    Flux Transfer Function
    # Add License
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

CC_A = CellData.Const.CC_A   # Current-collector area [m^2]
as = 3*Electrode.ϵ_s/Electrode.Rs # Specific interfacial surf. area
κ_eff = CellData.Const.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = CellData.Const.SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0 

#Beta's
β = @. Electrode.Rs*sqrt(s/Electrode.Ds)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max = Electrode.cs_max
cs0 = cs_max * θ
α = Electrode.α

#Current Flux Density
κ = Electrode.k_norm/Electrode.cs_max/ce0^(1-α)
j0 = κ*(ce0*(cs_max-cs0))^(1-α)*cs0^α

#Resistances
Rtot = R*CellData.Const.T/(j0*F^2) + Electrode.RFilm

#∂Uocp_Def
∂Uocp_elc = CellData.Const.∂Uocp(Def,θ)/cs_max

cse_res = -3/(as*F*Electrode.L*CC_A*Electrode.Rs) #Residual Variable for Pole Removal - eq. 4.45
ν = @. Electrode.L*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(β)/(tanh(β)-β)))) #Condensing Variable - eq. 4.13

cse_tf = @. (ν*Electrode.Rs*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))*tanh(β))/(as*F*Electrode.L*CC_A*Electrode.Ds*(κ_eff+σ_eff)*sinh(ν)*(tanh(β)-β)) #Transfer Function - eq. 4.17
cse_tf =  cse_tf.-cse_res./s  #Pole removal - eq. 4.43

zero_tf = @. (5*as*Electrode.Ds*F*Electrode.L^2*(κ_eff*(2-6*z+3*z^2)+(3*z^2-1)*σ_eff)-6*∂Uocp_elc*Electrode.Rs*κ_eff*σ_eff*κ_eff)/(30*CC_A*as*Electrode.Ds*∂Uocp_elc*F*Electrode.L*κ_eff*σ_eff) #For s = 0 / Wolfram Alpha
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
        println("as:C_se:Pos",as)
        println("Rtot:C_se:Pos",Rtot)
        println("σ_eff:C_se:Pos",σ_eff)
        println("∂Uocp_elc:C_se:Pos",∂Uocp_elc)
        println("Electrode.L:C_se:Pos",Electrode.L)
    end
end
cse_res = cse_res*ones(length(z))
return cse_tf, D, cse_res, D_term

end
