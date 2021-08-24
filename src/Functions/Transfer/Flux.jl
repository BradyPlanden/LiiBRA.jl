@inline function Flux(CellData,s,z,Def)
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
β = @. Electrode.Rs*sqrt(s/Electrode.Ds )

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max = Electrode.cs_max
cs0 = cs_max * θ
α = Electrode.α

#Current Flux Density
if CellData.Const.CellTyp == "Doyle_94"
    κ = Electrode.k_norm/Electrode.cs_max/ce0^(1-α)
    j0 = κ*(ce0*(cs_max-cs0))^(1-α)*cs0^α
else
    j0 = Electrode.k_norm*(ce0*(cs_max-cs0))^(1-α)*cs0^α
end

#Resistances
Rtot = R*CellData.Const.T /(j0*F^2) + Electrode.RFilm

#∂Uocp_Def
∂Uocp_elc = CellData.Const.∂Uocp(Def,θ)/cs_max

#Condensing Variable
ν = @. Electrode.L*sqrt((as/σ_eff+as/κ_eff)/(Rtot+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(β)/(tanh(β)-β))))
ν_∞ = @. Electrode.L*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))

#Transfer Function
j_tf = @. ν*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))/(as*F*Electrode.L*CC_A*(κ_eff+σ_eff)*sinh(ν))
D = @. ν_∞*(σ_eff*cosh(ν_∞*z)+κ_eff*cosh(ν_∞*(z-1)))/(as*F*Electrode.L*CC_A*(κ_eff+σ_eff)*sinh(ν_∞))
D_term = "@. $ν_∞*($σ_eff*cosh($ν_∞*$z)+$κ_eff*cosh($ν_∞*($z-1)))/($as*$F*$(Electrode.L)*$CC_A*($κ_eff+$σ_eff)*sinh($ν_∞))"
zero_tf =ones(2)*1/(CellData.Const.CC_A*as*F*Electrode.L)
j_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]
res0 = zeros(length(z))

if Def == "Pos"
    j_tf = -j_tf
    D = -D
    D_term = "@. -$ν_∞*($σ_eff*cosh($ν_∞*$z)+$κ_eff*cosh($ν_∞*($z-1)))/($as*$F*$(Electrode.L)*$CC_A*($κ_eff+$σ_eff)*sinh($ν_∞))"
    if Debug == 1
        println("D:Flux:Pos",D)
        println("z:Flux:Pos",z)
        println("ν_∞:Flux:Pos",ν_∞)
        println("j0:Flux:Pos",j0)
        println("Rtot:Flux:Pos",Rtot)
        println("σ_eff:Flux:Pos",σ_eff)
        println("ν_∞:Flux:Pos",ν_∞)
        println("L:Flux:Pos",Electrode.L)
    end
else 
    if Debug == 1  
        println("D:Flux:Neg",D)
        println("z:Flux:Neg",z)
        println("ν_∞:Flux:Neg",ν_∞)
        println("D:Flux:Neg",D)
    end
end

return j_tf, D, res0, D_term

end
