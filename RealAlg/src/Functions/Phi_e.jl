@inline function Phi_e(CellData::Cell,s,z)
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
Ds_Neg = CellData.Neg.Ds       # Solid diffusivity [m^2/s]
Ds_Pos = CellData.Pos.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Geo.CC_A   # Current-collector area [m^2]
De = CellData.Const.De # Electrolyte Diffusivity
κ_eff_Neg = CellData.Const.κ*ϵ1^CellData.Neg.κ_brug
κ_eff_Sep = CellData.Const.κ*ϵ2^CellData.Sep.κ_brug
κ_eff_Pos = CellData.Const.κ*ϵ3^CellData.Pos.κ_brug
σ_eff_Neg = CellData.Neg.σ*ϵ1^CellData.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = CellData.Pos.σ*ϵ3^CellData.Pos.σ_brug #Effective Conductivity Pos
dln = CellData.Const.dln  #Electrolyte activity coefficient term (Rod. 17)
κ_D_eff = (2*R*T/F)*κ_eff_Neg*(1-t_plus)*(1+dln) #Diffision Effective Electrolyte Conductivity

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
Rct_neg = R*T/(j0_neg*F)^2
Rtot_neg = Rct_neg + CellData.Neg.RFilm

Rct_pos = R*T/(j0_pos*F)^2
Rtot_pos = Rct_pos + CellData.Pos.RFilm

#OCP derivative
∂Uocp_pos = ∂Uocp("Pos",θ_pos)
∂Uocp_neg = ∂Uocp("Neg",θ_neg)

ν_neg = @. Lneg*sqrt((as_neg/σ_eff_Neg+as_neg/κ_eff_Neg)/(Rtot_neg.+∂Uocp_neg*(Rs_Neg/(F*Ds_Neg)).*(tanh.(βn)./(tanh.(βn)-βn)))) #Condensing Variable - eq. 4.13
ν_neg_∞ = @. Lneg*sqrt(as_neg*((1/κ_eff_Neg)+(1/σ_eff_Neg))/(Rtot_neg))
ν_pos = @. Lpos*sqrt((as_pos/σ_eff_Pos+as_pos/κ_eff_Pos)/(Rtot_pos.+∂Uocp_pos*(Rs_Pos/(F*Ds_Pos)).*(tanh.(βp)./(tanh.(βp)-βp)))) #Condensing Variable - eq. 4.13
ν_pos_∞ = @. Lpos*sqrt(as_pos*((1/κ_eff_Pos)+(1/σ_eff_Pos))/(Rtot_pos))


# if tf = Rod17 #Consider splitting into different function calls.

#    for i < length(z)
#    pt = z[i]

#       if pt <= Lneg
#          ϕ_tf = @. Lneg*(σ_eff*(1-cosh(ν_neg*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-pt))-pt*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
#          zero_tf = @. -(Lneg*pt^2)/(2*as*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
#          #D_term  =  @. Lneg*(σ_eff*(1-cosh(ν_∞*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (Lneg*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞)
#          ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

#       elseif pt <= Lneg + Lsep
#          ϕ_tf1 = @. -pt*Lsep/(as*κ_eff)
#          ϕ_tf2 = @. -(κ_D_eff/κ_eff_Pos*ce0)*(cS1*(e^ΛS1*pt-1)+cS2*(e^ΛS1*pt-1))
#          ϕ_n_tf = @. Lneg*(σ_eff*(1-cosh(ν_neg*Lneg)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-Lneg))-Lneg*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
#          ϕ_tf = @. ϕ_n_tf + ϕ_tf1 + ϕ_tf2 
#          zero_tf = @. -(Lneg*pt^2)/(2*as*κ_eff) - (κ_d*Lneg)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
#          #D_term  = @. Lsep*(σ_eff*(1-cosh(ν_∞*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (Lsep*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) - pt*Lsep/(as*κ_eff)
#          ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

#       else
#          ϕ_s_tf1 = @. -pt*Lsep/(as*κ_eff)
#          ϕ_s_tf2 = @. -(κ_D_eff/κ_eff_Pos*ce0)*(cS1*(e^ΛS1*pt-1)+cS2*(e^ΛS1*pt-1))
#          ϕ_n_tf = @. Lneg*(σ_eff*(1-cosh(ν_neg*Lneg)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (Lneg*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-Lneg))-Lneg*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
#          ϕ_p_tf1 = ((pt-1)*Lpos/κ_eff_Pos*CC_A)+(as_pos*F*Lpos^2/κ_eff_Pos) * (j_Pos_1*(e^(ΛP1*pt)+(1-pt)*ΛP1*e^ΛP1-e^ΛP1)/ΛP1^2 + j_Pos_2*(e^(-ΛP1*pt)+(pt-1)*ΛP1*e^-ΛP1-e^-ΛP1)/ΛP1^2 + j_Pos_1*(e^(ΛP1*pt)+(1-pt)*ΛP1*e^ΛP1-e^ΛP1)/ΛP1^2 + j_Pos_3*(e^(ΛP2*pt)+(1-pt)*ΛP2*e^ΛP2-e^ΛP2)/ΛP2^2 + j_Pos_4*(e^(-ΛP2*pt)+(pt-1)*ΛP2*e^-ΛP2-e^-ΛP2)/ΛP2^2)
#          ϕ_p_tf2 = -(κ_D_eff/κ_eff_Pos*ce0)*(cp1*(e^ΛP1*pt-e^Λp1)+cp2*(e^ΛP1*pt-e^Λp1)+cp3*(e^ΛP2*pt-e^Λp2)+cp4*(e^ΛP2*pt-e^Λp2))
#          ϕ_tf = @. ϕ_n_tf + ϕ_s_tf1 + ϕ_s_tf2 + ϕ_p_tf1 + ϕ_p_tf2

#          #D_term  = @.  
#    end

ϕ_tf = fill(0.0 + im,(length(z),length(s)))
#ϕ_tf = convert(Complex{Float64},ϕ_tf)

D_term = fill(0.0 + im,(length(z),length(s)))
#D_term = convert(Complex{Float64},ϕ_tf)

# else
   for i = 1:length(z)
      pt = z[i]
   
         if pt <= Lneg
            ϕ_tf[i,:] = @. (Lneg*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg*pt/Lneg)) - pt*ν_neg*sinh(ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg) + (Lneg*(cosh(ν_neg)-cosh(ν_neg*(Lneg-pt)/Lneg)/(CC_A*κ_eff_Neg*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg))) #Lee. Eqn. 4.22
            zero_tf = @. -(pt^2)/(2*CC_A*κ_eff_Neg*Lneg)
            D_term[i,:]  =  @. (Lneg*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg*pt/Lneg)) - pt*ν_neg*sinh(ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg) + (Lneg*(cosh(ν_neg)-cosh(ν_neg*(Lneg-pt)/Lneg)/(CC_A*κ_eff_Neg*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg)))
            ϕ_tf[findall(s.==0),:] .= zero_tf
   
         elseif pt <= Lneg + Lsep
            ϕ_tf[i,:] = @. (Lneg - pt)/(CC_A*κ_eff_Sep) + (Lneg*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2))-ν_neg)/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) #Lee. Eqn. 4.23
            zero_tf = @. (2*κ_eff_Neg*Lneg-κ_eff_Sep*Lneg-2*κ_eff_Neg*pt)/(2*CC_A*κ_eff_Neg*κ_eff_Sep)
            D_term[i,:] .= @. (Lneg - pt)/(CC_A*κ_eff_Sep) + (Lneg*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg_∞/2))-ν_neg_∞)/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg_∞)
            ϕ_tf[findall(s.==0),:] .= zero_tf
   
         else
            ϕ_tf[i,:] = @. (Lsep/(CC_A*κ_eff_Sep)) + Lneg*(((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2))-ν_neg)/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) - Lpos*(1+(σ_eff_Pos/κ_eff_Pos)*cosh(ν_pos))/(CC_A*(κ_eff_Pos+σ_eff_Neg)*sinh(ν_pos)*ν_pos)  #Lee. Eqn. 4.24
            zero_tf = @. -(κ_eff_Sep*κ_eff_Pos*Lneg*Lpos)/(2*CC_A*κ_eff_Neg*κ_eff_Sep*κ_eff_Pos*Lpos) + (κ_eff_Neg*(-2*κ_eff_Pos*Lpos*Lsep+κ_eff_Sep*(Lneg+Lsep-pt)*(Lneg+2*Lpos+Lsep-pt)))/(2*CC_A*κ_eff_Neg*κ_eff_Sep*κ_eff_Pos*Lpos)
            D_term[i,:]  .= @. (Lsep/(CC_A*κ_eff_Sep)) + Lneg*(((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg_∞/2))-ν_neg_∞)/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg_∞) - Lpos*(1+(σ_eff_Pos/κ_eff_Pos)*cosh(ν_pos_∞))/(CC_A*(κ_eff_Pos+σ_eff_Neg)*sinh(ν_pos_∞)*ν_pos_∞)
            ϕ_tf[findall(s.==0),:] .= zero_tf
      end
 end
#println("D_term",size(D_term))
return ϕ_tf, D_term

end
