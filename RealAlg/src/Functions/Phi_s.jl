function Phi_s(CellData::Cell,FCall::FCalls,s,z,Def)
    """ 
    Solid Potential Transfer Function
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

L = Electrode.L #Electrode Length
T = CellData.Const.T  # Temperature
as = 3*Electrode.ϵ_s/Electrode.Rs # Specific interfacial surf. area
κ_eff = FCall.Kap.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = CellData.Const.Init_SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0

#Beta's
β = @. Electrode.Rs*sqrt(s/Electrode.Ds)

#Prepare for j0
cs0 = Electrode.cs_max * θ

#Current Flux Density
κ = Electrode.k_norm/Electrode.cs_max/CellData.Const.ce0^(1-Electrode.α)
j0 = κ*(CellData.Const.ce0*(Electrode.cs_max-cs0))^(1-Electrode.α)*cs0^Electrode.α

#Resistances
Rtot = R*T/(j0*F^2) + Electrode.RFilm

#∂Uocp_Def = UOCP(θ_Def)
∂Uocp_elc = ∂Uocp(Def,θ)/Electrode.cs_max

ν = @. L*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(β)/(tanh(β)-β)))) #Condensing Variable - eq. 4.13
ν_∞ = @. L*sqrt((as*(1/κ_eff)+(1/σ_eff))/(Rtot))

ϕ_tf = @. -L*(κ_eff*(cosh(ν)-cosh(z-1)*ν))/(CellData.Geo.CC_A*σ_eff*(κ_eff+σ_eff)*ν*sinh(ν))-L*(σ_eff*(1-cosh(z*ν)+z*ν*sinh(ν)))/(CellData.Geo.CC_A*σ_eff*(κ_eff+σ_eff)*ν*sinh(ν)) #Transfer Function - eq. 4.19
D_term = @. -L*(κ_eff*(cosh(ν_∞)-cosh(z-1)*ν_∞))/(CellData.Geo.CC_A*σ_eff*(κ_eff+σ_eff)*ν_∞*sinh(ν_∞))-L*(σ_eff*(1-cosh(z*ν_∞)+z*ν_∞*sinh(ν_∞)))/(CellData.Geo.CC_A*σ_eff*(κ_eff+σ_eff)*ν_∞*sinh(ν_∞)) # Contribution to D as G->∞
D_term_check = @. ((-2+z)*z*L)/(2*CellData.Geo.CC_A*σ_eff)
zero_tf = @. L*(z-2)*z/(2*CellData.Geo.CC_A*σ_eff)
ϕ_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]
res0 = zeros(length(z))

if Def == "Pos" #Double check this implementation
   ϕ_tf = -ϕ_tf
   D_term = -D_term
   if Debug == 1
      println("D_term:Phi_s:Pos:",D_term)
      println("D_term_check:Phi_s:Pos:",D_term_check)
      #println("z:Phi_s:Pos:",z)
      println("ν_∞:Phi_s:Pos:",ν_∞)
      #println("ν:Phi_s:Pos:",ν)
      println("θ:Phi_s:Pos:",θ)
      println("κ:Phi_s:Pos:",κ)
      println("Rtot:Phi_s:Pos:",Rtot)
      println("j0:Phi_s:Pos:",j0)
      println("κ_eff:Phi_s:Pos:",κ_eff)
      println("σ_eff:Phi_s:Pos:",σ_eff)
   end
else
   if Debug == 1
      println("D_term:Phi_s:Neg:",D_term)
      println("D_term_check:Phi_s:Neg:",D_term_check)
      println("ν_∞:Phi_s:Neg:",ν_∞)
      println("θ:Phi_s:Neg:",θ)
      println("κ:Phi_s:Neg:",κ)
      println("Rtot:Phi_s:Neg:",Rtot)
      println("j0:Phi_s:Neg:",j0)
      println("κ_eff:Phi_s:Neg:",κ_eff)
      println("σ_eff:Phi_s:Neg:",σ_eff)
      println("zero_tf:Phi_s:Neg:",zero_tf)
      #println("z:Phi_s:Neg:",z)
      #println("ν:Phi_s:Neg:",ν)
   end
end
return ϕ_tf, D_term, res0

end
