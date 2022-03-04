@inline function C_se(Cell,s,z,Def,cse_tf,D,res0)
""" 
Concentration Solid-Electrolyte Transfer Function

C_se(Cell,s,z,Def)

"""

if Def == "Pos"   
    Electrode = Cell.Pos #Electrode Length
else
    Electrode = Cell.Neg #Electrode Length
end


κ_eff = Cell.Const.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

#Defining SOC
θ = Cell.Const.SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0 

#Prepare for j0
cs0 = Electrode.cs_max * θ

#Current Flux Density
if Cell.Const.CellTyp == "Doyle_94"
    κ = Electrode.k_norm/Electrode.cs_max/Cell.Const.ce0^(1-Electrode.α)
    j0 = κ*(Cell.Const.ce0*(Electrode.cs_max-cs0))^(1-Electrode.α)*cs0^Electrode.α
else
    j0 = Electrode.k_norm*(Cell.Const.ce0*cs0*(Electrode.cs_max-cs0))^(1-Electrode.α)
end

#Resistances
Rtot = R*Cell.Const.T/(j0*F^2) + Electrode.RFilm
#Rtot = R*Cell.Const.T/(j0*Cell.Const.CC_A*F) + Electrode.RFilm

#∂Uocp_Def
∂Uocp_elc = Cell.Const.∂Uocp(Def,θ)/Electrode.cs_max

res0 .= -3/(Electrode.as*F*Electrode.L*Cell.Const.CC_A*Electrode.Rs)*ones(length(z)) #Residual Variable for Pole Removal - eq. 4.45
ν = @. Electrode.L*sqrt((Electrode.as/σ_eff+Electrode.as/κ_eff)/(Rtot+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(Electrode.β)/(tanh(Electrode.β)-Electrode.β)))) #Condensing Variable - eq. 4.13

cse_tf .= @. ν*Electrode.Rs*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))*tanh(Electrode.β)/(Electrode.as*F*Electrode.L*Cell.Const.CC_A*Electrode.Ds*(κ_eff+σ_eff)*sinh(ν)*(tanh(Electrode.β)-Electrode.β))-res0/s #Transfer Function - eq. 4.17
#cse_tf .= @. cse_tf-res0t/s  #Pole removal - eq. 4.43
zero_tf = @. (5*Electrode.as*Electrode.Ds*F*Electrode.L^2*(κ_eff*(2-6*z+3*z^2)+(3*z^2-1)*σ_eff)-6*∂Uocp_elc*Electrode.Rs*κ_eff*σ_eff)/(30*Cell.Const.CC_A*Electrode.as*Electrode.Ds*∂Uocp_elc*F*Electrode.L*κ_eff*σ_eff) #For s = 0 / Wolfram Alpha
cse_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]
D .= zeros(length(z))

if Def == "Pos"
    cse_tf .= -cse_tf
    res0 .= -res0
end

if abs.(cse_tf[:,1]) > abs.(cse_tf[:,2])*10
    cse_tf[:,1] = cse_tf[:,2]*10
end

end
