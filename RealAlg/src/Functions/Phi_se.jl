function Phi_se(CellData::Cell,s,z,Def)
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
β = permutedims(β)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max = Electrode.cs_max
cs0 = cs_max * θ
α = Electrode.α


j0 = Electrode.k_norm*(ce0*(cs_max-cs0))^(1-α)*cs0^α #Exchange Current Density

#Resistances
Rct = R*T/(j0*F)^2
Rtot = Rct + Electrode.RFilm


∂Uocp_elc = ∂Uocp(Def,θ)/cs_max #Open Circuit Potential Partial

res0 = -3*(∂Uocp_elc)/(as*F*L*CC_A*Rs) # residual for pole removal

ν = @. L*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Rs/(F*Ds))*(tanh(β)/(tanh(β)-β)))) #Condensing Variable - eq. 4.13
ν_∞ = @. L*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))

println("ν:",size(ν))
println("res0:",res0)
println("∂Uocp_elc:",∂Uocp_elc)

ϕ_tf = @. L/(CC_A*ν*sinh(ν)*((1/κ_eff)*cosh(ν*z)+(1/σ_eff)*cosh(ν*(z-1))))-res0  #Transfer Function - eq. 4.14
zero_tf = @. (6*κ_eff*(5*Ds*F*Rtot-∂Uocp_elc*Rs)*σ_eff)/(30*CC_A*as*Ds*F*κ_eff*σ_eff) + 5*as*Ds*F*L^2*(σ_eff*(-1+3*z^2)+κ_eff*(2-6*z+3*z^2)/(30*CC_A*as*Ds*F*κ_eff*L*σ_eff))
D_term = @. L/(CC_A*ν_∞*sinh(ν_∞)*((1/κ_eff)*cosh(ν_∞*z)+(1/σ_eff)*cosh(ν_∞*(z-1)))) # Contribution to D as G->∞
ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

if Def == "Pos" #Double check this implementation
   ϕ_tf = -ϕ_tf
end

return ϕ_tf, D_term

end
