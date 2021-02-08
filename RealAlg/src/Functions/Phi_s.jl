function Phi_s(CellData::Cell,s,z,Def)
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
T = CellData.Const.T      # Temperature
t_plus = CellData.Const.t_plus  # Transference Number
ζ = (1-t_plus)/F    #Simplifying Variable
Rs = Electrode.Rs       # Particle radius [m]
Ds = Electrode.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Geo.CC_A   # Current-collector area [m^2]
as = 3*Electrode.ϵ_s/Rs # Specific interfacial surf. area
κ_eff = CellData.Const.κ*ϵ1^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*ϵ1^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = CellData.Const.Init_SOC * (Electrode.θ_max-Electrode.θ_min) + Electrode.θ_min 

#Beta's
β = @. Rs*sqrt(s/Ds)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max = Electrode.cs_max
cs0 = cs_max * θ
α = Electrode.α

#Current Flux Density
j0 = Electrode.k_norm*(ce0*(cs_max-cs0))^(1-α)*cs0^α

#Resistances
Rct = R*T/(j0*F)^2
Rtot = Rct + Electrode.RFilm

#∂Uocp_Def = UOCP(θ_Def)
∂Uocp_elc = ∂Uocp(Def,θ)

ν = @. L*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Rs/(F*Ds)).*(tanh.(β)./(tanh.(β)-β)))) #Condensing Variable - eq. 4.13
ν_∞ = @. L*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))

ϕ_tf = @. L*(κ_eff*(cosh(ν)-cosh(z'-1)*ν)/(as*σ_eff*(κ_eff+σ_eff)*ν*sinh(ν)-L*(σ_eff*(1-cosh(z'*ν)+z'*ν*sinh(ν)))/(as*σ_eff*(κ_eff+σ_eff)*ν*sinh(ν)))) #Transfer Function - eq. 4.19
D_term = @. L*(κ_eff*(cosh(ν_∞)-cosh(z'-1)*ν_∞)/(as*σ_eff*(κ_eff+σ_eff)*ν_∞*sinh(ν_∞)-L*(σ_eff*(1-cosh(z'*ν_∞)+z'*ν_∞*sinh(ν_∞)))/(as*σ_eff*(κ_eff+σ_eff)*ν_∞*sinh(ν_∞)))) # Contribution to D as G->∞


if Def == "Pos" #Double check this implementation
   ϕ_tf = -ϕ_tf
end

return ϕ_tf, D_term

end
