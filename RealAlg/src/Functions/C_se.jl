@inline function C_se(CellData::Cell,FCall::FCalls,s,z,Def)
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
σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

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
κ = Electrode.k_norm/Electrode.cs_max/ce0^(1-α)
j0 = κ*(ce0*(cs_max-cs0))^(1-α)*cs0^α

#Resistances
Rct = R*T/(j0*F)^2
Rtot = Rct + Electrode.RFilm

#∂Uocp_Def = UOCP(θ_Def)
∂Uocp_elc = ∂Uocp(Def,θ)

cse_res = -3/(as*F*L*CC_A*Ds) #Residual Variable for Pole Removal - eq. 4.45
ν = @. L*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Rs/(F*Ds)).*(tanh.(β)./(tanh.(β)-β)))) #Condensing Variable - eq. 4.13
cse_tf = @. ν*Rs*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))/(as*F*L*CC_A*Ds*(κ_eff+σ_eff)*sinh(ν)) #Transfer Function - eq. 4.17
cse_tf = @. cse_tf-cse_res #Pole removal - eq. 4.43
zero_tf = @. (5*as*Ds*F*L^2*((2-6z+3z^2)+(3z^2-1)*σ_eff)-6*∂Uocp_elc*Rs*κ_eff*σ_eff*κ_eff)/(30*CC_A*as*Ds*∂Uocp_elc*F*L*κ_eff*σ_eff) #For s = 0 / Wolfram Alpha
cse_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

if Def == "Pos" #Double check this implementation
    cse_tf = -cse_tf
end
D_term = zeros(length(z))
return cse_tf, D_term

end
