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


T = CellData.Const.T      # Temperature
t_plus = CellData.Const.t_plus  # Transference Number
ζ = (1-t_plus)/F    #Simplifying Variable
Rs = Electrode.Rs       # Particle radius [m]
Ds_Neg = CellData.Neg.Ds       # Solid diffusivity [m^2/s]
Ds_Pos = CellData.Pos.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Geo.CC_A   # Current-collector area [m^2]
as = 3*Electrode.ϵ_s/Rs # Specific interfacial surf. area
De = CellData.Const.De # Electrolyte Diffusivity
κ_eff_Neg = CellData.Const.κ*ϵ1^CellData.Neg.κ_brug
κ_eff_Pos = CellData.Const.κ*ϵ3^CellData.Pos.κ_brug
σ_eff_Neg = CellData.Neg.σ*ϵ1^CellData.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = CellData.Pos.σ*ϵ3^CellData.Pos.σ_brug #Effective Conductivity Pos
dln = 3  #Electrolyte activity coefficient term (Rod. 17)
κ_D_eff = (2*R*T/F)*κ_eff(1-t_plus)(1+dln) #Diffision Effective Electrolyte Conductivity

#Defining SOC
θ_neg = CellData.Const.Init_SOC * (CellData.Neg.θ_max-CellData.Neg.θ_min) + CellData.Neg.θ_min 
θ_pos = CellData.Const.Init_SOC * (CellData.Pos.θ_max-CellData.Pos.θ_min) + CellData.Pos.θ_min 

#Beta's
βn = Rs_Neg.*sqrt.(s./Ds_Neg)
βp = Rs_Pos.*sqrt.(s./Ds_Pos)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max_neg = CellData.Neg.cs_max
cs0_neg = cs_max_neg * θ_neg
cs_max_pos = CellData.Pos.cs_max
cs0_pos = cs_max_pos * θ_pos
α_neg = CellData.Neg.α
α_pos = CellData.Pos.α

#Current Flux Density
j0_neg = CellData.Neg.k_norm*(ce0*(cs_max_neg-cs0_neg))^(1-α_neg)*cs0_neg^α_neg
j0_pos = CellData.Pos.k_norm*(ce0*cs_max_pos*cs0_pos)^(1-α_pos)*cs0_pos^α_pos

#Resistances
Rct = R*T/(j0*F)^2
Rtot = Rct + Electrode.RFilm

#OCP derivative
∂Uocp_pos = ∂Uocp('Pos',θ_pos)
∂Uocp_neg = ∂Uocp('Neg',θ_neg)

ν_neg = @. Lneg*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Rs/(F*Ds)).*(tanh.(β)./(tanh.(β)-β)))) #Condensing Variable - eq. 4.13
ν_neg_∞ = @. Lneg*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))
ν_pos = @. Lpos*sqrt((as/σ_eff+as/κ_eff)/(Rtot.+∂Uocp_elc*(Rs/(F*Ds)).*(tanh.(β)./(tanh.(β)-β)))) #Condensing Variable - eq. 4.13
ν_pos_∞ = @. Lpos*sqrt(as*((1/κ_eff)+(1/σ_eff))/(Rtot))


if tf = Rod17 #Consider splitting into different function calls.

   for i < length(z)
   pt = z[i]

      if pt <= Lneg
         ϕ_tf = @. Lneg*(σ_eff*(1-cosh(ν_neg*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-pt))-pt*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
         zero_tf = @. -(Lneg*pt^2)/(2*as*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
         D_term  =  @. Lneg*(σ_eff*(1-cosh(ν_∞*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (Lneg*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞)
         ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

      elseif pt <= Lneg + Lsep
         ϕ_tf1 = @. -pt*Lsep/(as*κ_eff)
         ϕ_tf2 = @. -(κ_D_eff/κ_eff_Pos*ce0)*(cS1*(e^ΛS1*pt-1)+cS2*(e^ΛS1*pt-1))
         ϕ_tf0 = @. Lneg*(σ_eff*(1-cosh(ν_neg*Lneg)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-Lneg))-Lneg*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
         ϕ_tf = @. ϕ_tf0 + ϕ_tf1 + ϕ_tf2 
         zero_tf = @. -(Lneg*pt^2)/(2*as*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
         D_term  = @. Lsep*(σ_eff*(1-cosh(ν_∞*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (Lsep*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) - pt*Lsep/(as*κ_eff)
         ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

      else
         ϕ_tf1 = ((pt-1)*Lpos/κ_eff_Pos*as)+()
         ϕ_tf2 = -(κ_D_eff/κ_eff_Pos*ce0)*(cp1*(e^ΛP1*pt-e^Λp1)+cp2*(e^ΛP1*pt-e^Λp1)+cp3*(e^ΛP2*pt-e^Λp2)+cp4*(e^ΛP2*pt-e^Λp2))
         ϕ_tf = @. σ_eff*sinh(z'*ν/Lneg)+κ_eff*sinh((Lneg-z')*ν/Lneg) / (as*(κ_eff+σ_eff)*sinh(ν)) + κ_eff/(as*((κ_eff+σ_eff)))
         D_term  = @.  
   end

else
   for i < length(z)
      pt = z[i]
   
         if pt <= Lneg
            ϕ_tf = @. (Lneg*(σ_eff/κ_eff)(1-cosh(ν_neg*pt/Lneg)) - pt*ν_neg*sinh(ν_neg))/(as*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(cosh(ν_neg)-cosh(ν_neg*(Lneg-pt)/Lneg)/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) #Lee. Eqn. 4.22
            zero_tf = @. -(Lneg*pt^2)/(2*as*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
            D_term  =  @. (Lneg*(σ_eff/κ_eff)(1-cosh(ν_neg*pt/Lneg)) - pt*ν_neg*sinh(ν_neg))/(as*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(cosh(ν_neg)-cosh(ν_neg*(Lneg-pt)/Lneg)/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
            ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]
   
         elseif pt <= Lneg + Lsep
            ϕ_tf = @. (Lneg - pt)/(as*κ_eff_sep) + (Lneg*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2))-ν_neg)/(as*(κ_eff_Neg+σ_eff_Neg)*ν_neg) #Lee. Eqn. 4.23
            zero_tf = @. -(Lneg*pt^2)/(2*as*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
            D_term  =@. (Lneg - pt)/(as*κ_eff_sep) + (Lneg*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg_∞/2))-ν_neg_∞)/(as*(κ_eff_Neg+σ_eff_Neg)*ν_neg_∞)
            ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]
   
         else
            ϕ_tf = @. (Lsep/(as*κ_eff_sep)) + Lneg*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2))-ν_neg)/(as*(κ_eff_Neg+σ_eff_Neg)*ν_neg) - Lpos*(1+(σ_eff_Pos/κ_eff_Pos)*cosh(ν_pos))/(as*(κ_eff_Pos+σ_eff_Neg)*sinh(ν_pos)*ν_pos))  #Lee. Eqn. 4.24
            D_term  = @.  
      end
end

if Def == "Pos" #Double check this implementation
  ϕ_tf = -ϕ_tf
end

return ϕ_tf, D_term

end
