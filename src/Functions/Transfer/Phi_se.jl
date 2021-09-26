@inline function Phi_se(CellData,s,z,Def)
    """ 
    Solid-Electrolyte Potential Transfer Function
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

ζ = (1-CellData.Const.t_plus)/F    #Simplifying Variable
Rs = Electrode.Rs       # Particle radius [m]
Ds = Electrode.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Const.CC_A   # Current-collector area [m^2]
as = 3*Electrode.ϵ_s/Rs # Specific interfacial surf. area
κ_eff = CellData.Const.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = CellData.Const.SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0

#Beta's
β = @. Rs*sqrt(s/Ds)

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
   j0 = Electrode.k_norm*(ce0*(cs0/cs_max*(1-cs0/cs_max)))^(1-α)
end

#Resistance
Rtot = R*CellData.Const.T/(j0*F^2) + Electrode.RFilm

∂Uocp_elc = CellData.Const.∂Uocp(Def,θ)/cs_max #Open Circuit Potential Partial
res0 = @. -3*∂Uocp_elc/(as*F*Electrode.L*CC_A*Rs) # residual for pole removal

ν = @. Electrode.L*sqrt((as/σ_eff+as/κ_eff)/(Rtot+∂Uocp_elc*(Rs/(F*Ds))*(tanh(β)/(tanh(β)-β)))) #Condensing Variable - eq. 4.13
ν_∞ = @. Electrode.L*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))

ϕ_tf = @. Electrode.L/(CC_A*ν*sinh(ν))*((1/κ_eff)*cosh(ν*z)+(1/σ_eff)*cosh(ν*(z-1))) #Transfer Function - eq. 4.14
ϕ_tf = @. ϕ_tf - res0./s

zero_tf = @. (6*(5*Ds*F*Rtot-∂Uocp_elc*Rs)*σ_eff)/(30*CC_A*as*Ds*F*σ_eff*Electrode.L) + (5*as*Ds*F*Electrode.L^2*(σ_eff*(-1+3*z^2)+κ_eff*(2-6*z+3*z^2)))/(30*CC_A*as*Ds*F*σ_eff*κ_eff*Electrode.L)
D = @. Electrode.L/(CC_A*ν_∞*sinh(ν_∞))*((1/κ_eff)*cosh(ν_∞*z)+(1/σ_eff)*cosh(ν_∞*(z-1))) # Contribution to D as G->∞
D_term = "@. $(Electrode.L)/($CC_A*$ν_∞*sinh($ν_∞))*((1/$κ_eff)*cosh($ν_∞*$z)+(1/$σ_eff)*cosh($ν_∞*($z-1)))"
ϕ_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]

if Def == "Pos" #Double check this implementation
   ϕ_tf = -ϕ_tf
   D = -D
   D_term = "@. -$(Electrode.L)/($CC_A*$ν_∞*sinh($ν_∞))*((1/$κ_eff)*cosh($ν_∞*$z)+(1/$σ_eff)*cosh($ν_∞*($z-1)))"
   if Debug == 1
      println("ϕ_tf:Phi_se:Pos",ϕ_tf[:,1])
      println("D:Phi_se:Pos",D)
      println("z:Phi_se:Pos",z)
      println("zero_tf:Phi_se:Pos",zero_tf)
      println("ν_∞:Phi_se:Pos",ν_∞)
      println("j0:Phi_se:Pos",j0)
      println("Rtot:Phi_se:Pos",Rtot)
      println("σ_eff:Phi_se:Pos",σ_eff)
      #println("res0:Phi_se:Pos",res0[:,1])
      println("∂Uocp_elc:Phi_se:Pos",∂Uocp_elc)
      println("L:Phi_se:Pos",Electrode.L)
  end
end
res0 = zeros(length(z))
return ϕ_tf, D, res0, D_term

end
