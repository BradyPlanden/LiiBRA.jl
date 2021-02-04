function Phi_e(CellData::Cell,s,z,Def)
   """ 
   Electrolyte Potential Transfer Function
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
As = 3*Electrode.ϵ_s/Rs # Specific interfacial surf. area
De = CellData.Const.De # Electrolyte Diffusivity
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

ν = @. L*sqrt((As/σ_eff+As/κ_eff)/(Rtot.+∂Uocp_elc*(Rs/(F*Ds)).*(tanh.(β)./(tanh.(β)-β)))) #Condensing Variable - eq. 4.13
ν_∞ = @. L*sqrt(As*((1/κ_eff)+(1/σ_eff))/(Rtot))


for i < length(z)
   pt = z[i]

      if pt <= Lneg
         ϕ_tf = @. Lneg*(σ_eff*(1-cosh(ν_neg*pt)/(As*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-pt))-pt*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞)
         zero_tf = @. -(Lneg*pt^2)/(2*As*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
         D_term  =  @. Lneg*(σ_eff*(1-cosh(ν_∞*pt)/(As*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (Lneg*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞)
         ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

      elseif pt <= Lpos && pt >= Lneg
         ϕ_tf = @. Lsep*(σ_eff*(1-cosh(ν_neg)/(As*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lsep*(κ_eff*cosh(ν_neg)-cosh(ν_neg)-sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) - pt*Lsep/(As*κ_eff)
         zero_tf = @. -(Lneg*pt^2)/(2*As*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
         D_term  = @. Lsep*(σ_eff*(1-cosh(ν_∞*pt)/(As*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (Lsep*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) - pt*Lsep/(As*κ_eff)
         ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]
      else
         ϕ_tf = @. σ_eff*sinh(z'*ν/Lneg)+κ_eff*sinh((Lneg-z')*ν/Lneg) / (As*(κ_eff+σ_eff)*sinh(ν)) + κ_eff/(As*((κ_eff+σ_eff)))
         D_term  = @.  
end

if Def == "Pos" #Double check this implementation
  ϕ_tf = -ϕ_tf
end

return ϕ_tf, D_term

end
