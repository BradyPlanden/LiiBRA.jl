@inline function j(CellData::Cell,FCall::FCalls,s,z,Def)
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

L = Electrode.L #Electrode Length
T = CellData.Const.T      # Temperature
t_plus = CellData.Const.t_plus  # Transference Number
ζ = (1-t_plus)/F    #Simplifying Variable
Rs = Electrode.Rs       # Particle radius [m]
Ds = Electrode.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Geo.CC_A   # Current-collector area [m^2]
as = 3*Electrode.ϵ_s/Rs # Specific interfacial surf. area
κ_eff = FCall.Kap.κ*ϵ1^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*ϵ1^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = CellData.Const.Init_SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0 

#Beta's
β = @. Rs*sqrt(s/Ds)
β = permutedims(β)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max = Electrode.cs_max
cs0 = cs_max * θ
α = Electrode.α

#Current Flux Density
j0 = Electrode.k_norm*(ce0*(cs_max-cs0))^(1-α)*cs0^α

#Resistances
Rct = R*T/(j0*F^2)
Rtot = Rct + Electrode.RFilm

#∂Uocp_Def = UOCP(θ_Def)
∂Uocp_elc = ∂Uocp(Def,θ)

#Condensing Variable
ν = @. L*sqrt((as/σ_eff+as/κ_eff)/(Rtot+∂Uocp_elc*(Rs/(F*Ds))*(tanh(β)/(tanh(β)-β))))
ν_∞ = @. L*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))

#Transfer Function
j_tf = @. ν*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))/(as*F*L*CC_A*(κ_eff+σ_eff)*sinh(ν))
D_term = @. ν_∞*(σ_eff*cosh(ν_∞*z)+κ_eff*cosh(ν_∞*(z-1)))/(as*F*L*CC_A*(κ_eff+σ_eff)*sinh(ν_∞))

if Def == "Pos" #Double check this implementation
    j_tf = -j_tf
    D_term = -D_term
end
#println("D_term:",D_term)
println("z:",z)
println("ν:",ν)
println("ν_∞:",ν_∞)
println("D_term:",D_term)

return j_tf, D_term

end
